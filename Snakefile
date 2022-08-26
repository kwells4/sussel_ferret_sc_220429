""" Snake pipeline for running cellranger with CITE-seq data """

# This has been borrowed and modified from pipelines written by Ryan Sheridan

# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
import subprocess
import glob
import os 
import re
from collections import defaultdict



# Parameters from config.yaml
RAW_DATA       = config["RAW_DATA"]
SAMPLES        = config["SAMPLES"]
RNA_SAMPLES    = config["RNA_SAMPLES"]
ADT_SAMPLES    = config["ADT_SAMPLES"]
VDJ_SAMPLES    = config["VDJ_SAMPLES"]
RESULTS        = config["RESULTS"]
GENOME         = config["GENOME"]
ADT_REF        = config["ADT_REF"]
VDJ_REF        = config["VDJ_REF"]
MAX_JOBS       = config["MAX_JOBS"]
LSF_TEMPLATE   = config["LSF_TEMPLATE"]
AGGR_GROUP     = config["AGGR_SAMPLES"]
CHEMISTRY      = config["CHEMISTRY"]
VELOCYTO_GROUP = config["VELOCYTO_GROUP"]
SANDBOX        = config["SANDBOX"]
STRINGTIE_GTF  = config["STRINGTIE_GTF"]

# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")



# Set sample/group names
RNA_SAMPLES     = [x.strip() for x in RNA_SAMPLES]
FASTQ_INFO      = "_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz"
FASTQ_SHORT     = "_S[0-9]+_R[12]_[0-9]+\.fastq\.gz"

SAMPLE_DICT_RNA = {SAMPLES[i]: RNA_SAMPLES[i] for i in range(len(SAMPLES))}

if ADT_SAMPLES:
    ADT_SAMPLES = [re.sub(", ", ",", x.strip()) for x in ADT_SAMPLES]
    #SAMPLES     = [x + "-" + re.sub(",", "_", y) for x, y in zip(RNA_SAMPLES, ADT_SAMPLES)]
    SAMPLE_DICT_ADT = {SAMPLES[i]: ADT_SAMPLES[i] for i in range(len(SAMPLES))}
    ADT_REF     = _check_path(ADT_REF)
else:
    SAMPLE_DICT_ADT = {SAMPLES[i]: "" for i in range(len(SAMPLES))}

if VDJ_SAMPLES:
    VDJ_SAMPLES = [x.strip() for x in VDJ_SAMPLES]
    SAMPLE_DICT_VDJ = {SAMPLES[i]: VDJ_SAMPLES[i] for i in range(len(SAMPLES))}
    VDJ_REF     = _check_path(VDJ_REF)
else:
    SAMPLE_DICT_VDJ = {SAMPLES[i]: "" for i in range(len(SAMPLES))}

print(CHEMISTRY)

# Make a chemistry dictionary
for i in range(len(SAMPLES)):
    if(CHEMISTRY is not None):
        if not CHEMISTRY.get(SAMPLES[i]):
            CHEMISTRY[SAMPLES[i]] = "auto"
    else:
        CHEMISTRY = {}
        CHEMISTRY[SAMPLES[i]] = "auto"

print(CHEMISTRY)

# Check/set directory/file paths
if len(RAW_DATA) == len(SAMPLES):
    RAW_DATA_DICT = {SAMPLES[i]: RAW_DATA[i] for i in range(len(SAMPLES))}
elif len(RAW_DATA) == 1:
    RAW_DATA_DICT = {SAMPLES[i]: RAW_DATA[0] for i in range(len(SAMPLES))}
else:
    sys.exit("RAW_DATA must either be the same length of samples or length 1")

if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)
RESULTS = _check_path(RESULTS)
GENOME  = _check_path(GENOME)

FASTQ_DIR = RESULTS + "/fastqs"
if not os.path.exists(FASTQ_DIR):
    os.makedirs(FASTQ_DIR)

if LSF_TEMPLATE:
    LSF_TEMPLATE = _check_path(LSF_TEMPLATE)
else:
    LSF_TEMPLATE = "lsf"

if not AGGR_GROUP:
    AGGR_GROUP = "none"

if not VELOCYTO_GROUP:
    VELOCYTO_GROUP = "none"

# Final output files
rule all:
    input:
        expand(
            "{results}/logs/{sample}_csv_done.out",
            results = RESULTS, sample = SAMPLES
        ),
        expand(
            "{results}/logs/{sample}_count_done.out",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
            "{results}/logs/{group}_csv_aggr_done.out",
            results = RESULTS, group = AGGR_GROUP
            ),
        expand(
            "{results}/logs/{group}_cellranger_aggr_done.out",
            results = RESULTS, group = AGGR_GROUP
            ),
        # Make bigwigs
        expand(
            "{results}/bigwig_trim/transfer_done.txt",
            results = RESULTS
            ),

        # Run stringtie
        expand(
            "{results}/stringtie/stringtie_combined.gtf",
            results = RESULTS
            )
        # expand(
        #     "{results}/logs/{sample}_velocyto_done.out",
        #     results = RESULTS, sample = SAMPLES
        # ),
        # expand(
        #     "{results}/logs/{group}_velocyto_combined.out",
        #     results = RESULTS, group = VELOCYTO_GROUP
        #     )

include: "src/rules/cellranger_multi.snake"
include: "src/rules/velocyto.snake"
include: "src/rules/make_bigwig.snake"
include: "src/rules/run_stringtie.snake"
