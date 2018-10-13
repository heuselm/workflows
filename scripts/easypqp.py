import sys
import os
import re
import operator
import pandas as pd

# alignment
from sklearn import preprocessing
import statsmodels.api as sm
from scipy.interpolate import interp1d

# mzML parsing
import pyopenms as po

def parse_tidemods(mods, term):
  mdb = po.ModificationsDB()
  modification_dictionary = []
  if mods != "":
    for mod in mods.split(','):
      if '+' in mod:
        sign = '+'
      elif '-' in mod:
        sign = '-'
      else:
        sys.exit("Error: Could not parse modifications.")

      # Check specificity
      specificity, mass = mod.split(sign)
      # Round mass
      mass = round(float(mass),2)
      # Treat static modifications differently
      staticmod = True
      if re.findall('\d+', specificity) != []:
        specificity = ''.join([i for i in specificity if not i.isdigit()])
        staticmod = False

      # Replace wild card
      if specificity == 'X':
        specificity = 'ARNDCEQGHILKMFPSTWYV'

      # Iterate over each residue
      for residue in specificity:
        if staticmod:
          residue_tide = residue
        else:
          if sign == '+':
            residue_tide = residue + "[" + str(mass) + "]"
          else:
            residue_tide = residue + "[" + sign + str(mass) + "]"

        if term == 'c':
          rm = mdb.getBestModificationByDiffMonoMass(mass, 0.02, residue, po.ResidueModification.TermSpecificity.C_TERM)
        elif term == 'n':
          rm = mdb.getBestModificationByDiffMonoMass(mass, 0.02, residue, po.ResidueModification.TermSpecificity.N_TERM)
        elif term == '':
          rm = mdb.getBestModificationByDiffMonoMass(mass, 0.02, residue, po.ResidueModification.TermSpecificity.ANYWHERE)

        if term == 'c':
          residue_unimod = residue + ".(UniMod:" + str(rm.getUniModRecordId()) + ")"
        elif term == 'n':
          residue_unimod = ".(UniMod:" + str(rm.getUniModRecordId()) + ")" + residue
        elif term == '':
          residue_unimod = residue + "(UniMod:" + str(rm.getUniModRecordId()) + ")"

        modification_dictionary.append({'tide': residue_tide, 'unimod': residue_unimod, 'unmodified': residue})
    return pd.DataFrame(modification_dictionary)
  else:
    return pd.DataFrame()

def remove_brackets(x):
  return re.sub(r'\[.*\]', '', x)

def read_percolator_psms(infile, runindex, fdr_threshold):
  table = pd.read_csv(infile, sep="\t")
  runs = pd.read_csv(runindex, sep=" ")

  psms = pd.merge(runs, table, on='file_idx')

  # Select proteotypic peptides only
  psms = psms[~psms['protein id'].str.contains(',')]

  # Select psms below q-value threshold
  psms = psms[psms['percolator q-value'] < fdr_threshold]

  return psms

def read_percolator_peptides(infile, fdr_threshold):
  peptides = pd.read_csv(infile, sep="\t")

  # Select peptides below q-value threshold
  peptides = peptides[peptides["percolator q-value"] < fdr_threshold]['sequence'].unique()

  return peptides

def read_percolator_proteins(infile, fdr_threshold):
  proteins = pd.read_csv(infile, sep="\t")

  # Select proteins below q-value threshold
  proteins = proteins[proteins["q-value"] < fdr_threshold]['ProteinId'].unique()

  return proteins

def lowess(run, reference_run):
  dfm = pd.merge(run, reference_run[['sequence','charge','irt']])

  print("INFO: Peptide overlap between run and reference: %s." % dfm.shape[0])

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm['irt'], dfm['rt'], frac=.66)
  lwf_x = list(zip(*lwf))[0]
  lwf_y = list(zip(*lwf))[1]
  lwi = interp1d(lwf_x, lwf_y, bounds_error=False)

  # Apply lowess model
  run['irt'] = lwi(run['rt'])

  return run

def read_mzxml(mzxml_path, scan_ids):
  fh = po.MzXMLFile()
  fh.setLogType(po.LogType.CMD)
  input_map = po.MSExperiment()
  fh.load(mzxml_path, input_map)

  rt_list = []
  peaks_list = []
  for scan_id in scan_ids:

    spectrum = input_map.getSpectrum(int(scan_id) - 1)
    rt_list.append(spectrum.getRT())

    product_mzs = []
    intensities = []
    for peak in spectrum:
      product_mzs.append(peak.getMZ())
      intensities.append(peak.getIntensity())

    peaks = pd.DataFrame({'filename': mzxml_path,'product_mz': product_mzs, 'intensity': intensities})
    peaks['precursor_mz'] = spectrum.getPrecursors()[0].getMZ()
    peaks['scan'] = scan_id
    peaks_list.append(peaks)

  transitions = pd.concat(peaks_list)
  rts = pd.DataFrame({'filename': mzxml_path, 'scan': scan_ids, 'rt': rt_list})
  return transitions, rts

# Parse input arguments
pepxmls = []
mzxmls = []
tidemods = sys.argv[1]
tidemods_cterm = sys.argv[2]
tidemods_nterm = sys.argv[3]
percolator_runindex = sys.argv[4]
percolator_psm_report = sys.argv[5]
psm_fdr_threshold = float(sys.argv[6])
percolator_peptide_report = sys.argv[7]
peptide_fdr_threshold = float(sys.argv[8])
percolator_protein_report = sys.argv[9]
protein_fdr_threshold = float(sys.argv[10])
outputdir = sys.argv[11] + "/"

for arg in sys.argv[12:]:
  if 'mzXML' in arg:
    mzxmls.append(arg)

# Parse Percolator reports
peptides = read_percolator_peptides(percolator_peptide_report, peptide_fdr_threshold)
print("INFO: Unique peptides below q-value threshold (%s): %s" % (peptide_fdr_threshold, len(peptides)))
proteins = read_percolator_proteins(percolator_protein_report, protein_fdr_threshold)
print("INFO: Unique proteins below q-value threshold (%s): %s" % (protein_fdr_threshold, len(proteins)))

psms = read_percolator_psms(percolator_psm_report, percolator_runindex, psm_fdr_threshold)

# Filter psms to all thresholds
psms = psms[(psms['sequence'].isin(peptides)) & (psms['protein id'].isin(proteins))]
print("INFO: Filtered redundant psms below q-value threshold (%s): %s" % (psm_fdr_threshold, psms.shape[0]))

# Patch Tide modifications
mods_nterm = parse_tidemods(tidemods_nterm,'n')
mods_cterm = parse_tidemods(tidemods_cterm,'c')
mods_residue = parse_tidemods(tidemods,'')

mods = pd.concat([mods_nterm, mods_cterm, mods_residue])

psms['unmodified_sequence'] = psms['sequence']
for idx, modification in mods.iterrows():
  print("Replace Tide modification '%s' with UniMod modification '%s'" % (modification['tide'], modification['unimod']))
  psms['unmodified_sequence'] = psms['unmodified_sequence'].str.replace(re.escape(modification['tide']), modification['unmodified'])
  psms['sequence'] = psms['sequence'].str.replace(re.escape(modification['tide']), modification['unimod'])

# Parse mzXML to retrieve peaks and metadata and store results in peak files
transitions = {}
rts = []
for mzxml in mzxmls:
  scans = psms[psms['filename'] == mzxml]['scan']
  print("INFO: Parsing file %s and extracting %s psms." % (mzxml, scans.shape[0]))
  transitions, rtdf = read_mzxml(mzxml, scans)
  transitions.to_pickle(outputdir + os.path.basename(os.path.splitext(mzxml)[0]+"_peaks.pkl"))
  rts.append(rtdf)
rtr = pd.concat(rts)

psms = pd.merge(psms, rtr, on=['filename','scan'])

# Generate set of best replicate identifications per run
psmsr = psms.loc[psms.groupby(['filename','sequence','charge'])['percolator q-value'].idxmin()].sort_index()

# Select reference run
psmsr_stats = psmsr.groupby('filename')[['sequence']].count().reset_index()
print(psmsr_stats)
reference_run_filename = psmsr_stats.loc[psmsr_stats['sequence'].idxmax()]['filename']

reference_run = psmsr[psmsr['filename'] == reference_run_filename]
align_runs = psmsr[psmsr['filename'] != reference_run_filename]

# Normalize RT of reference run
min_max_scaler = preprocessing.MinMaxScaler()
reference_run['irt'] = min_max_scaler.fit_transform(reference_run[['rt']])*100

# Normalize RT of all runs against reference
aligned_runs = align_runs.groupby('filename').apply(lambda x: lowess(x, reference_run)).dropna()
psmsa = pd.concat([reference_run, aligned_runs]).reset_index(drop=True)

# Generate set of non-redundant global best replicate identifications
psmsb = psmsa.loc[psmsa.groupby(['sequence','charge'])['percolator q-value'].idxmin()].sort_index()

# Write peak files
for mzxml in mzxmls:
  print("INFO: Generating peak reports for file %s." % mzxml)
  meta_global = psmsb[psmsb['filename'] == mzxml]

  transitions = pd.read_pickle(outputdir + os.path.basename(os.path.splitext(mzxml)[0]+"_peaks.pkl"))
  os.remove(outputdir + os.path.basename(os.path.splitext(mzxml)[0]+"_peaks.pkl"))

    # Generate run-specific PQP files for OpenSWATH alignment
  if "_Q1" in mzxml:
    meta_run = psmsa[psmsa['filename'] == mzxml]
    run_pqp = pd.merge(meta_run, transitions, on='scan')[['precursor_mz','product_mz','intensity','irt','protein id','unmodified_sequence','sequence','charge']]
    run_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
    run_pqp_path = outputdir + os.path.basename(os.path.splitext(mzxml)[0]+"_run_peaks.tsv")
    run_pqp.to_csv(run_pqp_path, sep="\t", index=False)

  # Generate global non-redundant PQP files
  global_pqp = pd.merge(meta_global, transitions, on='scan')[['precursor_mz','product_mz','intensity','irt','protein id','unmodified_sequence','sequence','charge']]
  global_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
  global_pqp_path = outputdir + os.path.basename(os.path.splitext(mzxml)[0]+"_global_peaks.tsv")
  global_pqp.to_csv(global_pqp_path, sep="\t", index=False)
