#!/bin/bash

## Directories
# HotSpot output directory (use absolute/complete directory location)
HD=$(pwd)/hotspot
# BAM file directory
BD=bam_3.2M
# additional BAM file suffix (use .rmd.srt if the files are .rmd.srt.bam, for example)
SUFFIX=.rmd.srt

# Make hotspot directories
mkdir $HD
mkdir ${HD}/tokens
mkdir ${HD}/output
mkdir ${HD}/peaks
mkdir ${HD}/regions

# HotSpot distr location
HSD=/data/mct-t200/T502/Lab_Notebook/tools/hotspot-distr

# Check genome and read length compatibility at
# http://www.uwencode.org/proj/hotspot/

# Genome
GENOME=mm10
# Read Length
RL=50

# Retrieve chromInfo and Mappability files from the UW ENCODE project web server
# If they aren't already present
CF=chromInfo.${GENOME}.bed
MF=${GENOME}.K${RL}.mappable_only.bed

if [ ! -f ${HD}/${CF} ]; then
    wget -O ${HD}/${CF} http://www.uwencode.org/proj/hotspot/${CF}
fi

if [ ! -f ${HD}/${MF} ]; then
    wget -O ${HD}/${MF} http://www.uwencode.org/proj/hotspot/${MF}
fi

# Build Hotspot Tokens
# See HotSpot's hotspot-distr/pipeline-scripts/test/runall.tokens.txt
# for a description of these settings
# cat EOM must not have excess white space.

for F in ${BD}/*.bam
do
NODIR=$(basename ${F})
BASE=`echo "$NODIR" | cut -d'.' -f1`

FILE=${HD}/tokens/${BASE}.txt
echo $FILE

cat <<EOM >$FILE
[script-tokenizer]
_TAGS_ = ${BD}/${BASE}${SUFFIX}.bam
_USE_INPUT_ = F
_INPUT_TAGS_ = 

_GENOME_ = mm10
_K_ = 50
_CHROM_FILE_ = ${HD}/${CF}
_MAPPABLE_FILE_ = ${HD}/${MF}

_DUPOK_ = T
_FDRS_ = "0.01"
_DENS_:

_OUTDIR_ = ${HD}/output/${BASE}
_RANDIR_ = ${HD}/output/${BASE}

_OMIT_REGIONS_:

_CHECK_ = T
_CHKCHR_ = chrX

_HOTSPOT_ = ${HSD}/hotspot-deploy/bin/hotspot
_CLEAN_ = T

_PKFIND_BIN_ = ${HSD}/hotspot-deploy/bin/wavePeaks
_PKFIND_SMTH_LVL_ = 3

_SEED_ = 101
_THRESH_ = 2
_WIN_MIN_ = 200
_WIN_MAX_ = 300
_WIN_INCR_ = 50
_BACKGRD_WIN_ = 50000
_MERGE_DIST_ = 150
_MINSIZE_ = 10
EOM
done


# Run Hotspot 
## NOTE the following section of this script is adapted from the runhotspot script 
## provided with HotSpot, and was not written by the authors of Gray, et al.
## It is included here for better integration with this analysis pipeline through
## use of shared directory variables.
#
#
## The folowing citation information for this code is from the HotSpot README file:
#
#
# Hotspot was originally conceived and written by Mike Hawrylycz, and is
# now maintained by Bob Thurman, University of Washington, with
# contributions by Eric Haugen, University of Washington, and especially
# Scott Kuehn.  Although there is no stand-alone publication for hotspot,
# the algorithm is described in detail in 
#
# Sam John et al., Chromatin accessibility pre-determines glucocorticoid
# receptor binding patterns, Nature Genetics 43, 264-268 
#
# The above should therefore serve as the primary citation for hotspot.
#
# This distribution is available via the uwencode website, at
#
#     http://uwencode.org/software/hotspot
#
# Bob Thurman
# rthurman@uw.edu
# 27 Jun 2013


set -e -o pipefail

scriptTokBin=${HSD}/ScriptTokenizer/src/script-tokenizer.py
pipeDir=${HSD}/pipeline-scripts

for tokenFile in ${HD}/tokens/*.txt
do

    ## Do everything, including badspots and final cleanup
    scripts="$pipeDir/run_badspot
        $pipeDir/run_make_lib
        $pipeDir/run_wavelet_peak_finding
        $pipeDir/run_10kb_counts
        $pipeDir/run_generate_random_lib
        $pipeDir/run_pass1_hotspot
        $pipeDir/run_pass1_merge_and_thresh_hotspots
        $pipeDir/run_pass2_hotspot
        $pipeDir/run_rescore_hotspot_passes
        $pipeDir/run_spot
        $pipeDir/run_thresh_hot.R
        $pipeDir/run_both-passes_merge_and_thresh_hotspots
        $pipeDir/run_add_peaks_per_hotspot
        $pipeDir/run_final"

    $scriptTokBin \
        --clobber \
        --output-dir=${HD}/output/ \
        $tokenFile \
        $scripts

    for script in $scripts
    do
        ${HD}/output/$(basename $script).tok
    done

    # Consolidate results
    NODIR=$(basename ${tokenFile})
    BASE=`echo "$NODIR" | cut -d'.' -f1`
    mv ./${BASE}*final ${HD}/output/${BASE}_final
    cp ${HD}/output/${BASE}_final/${BASE}${SUFFIX}.fdr0.01.pks.bed ${HD}/peaks/${BASE}.narrowPeak
    cp ${HD}/output/${BASE}_final/${BASE}${SUFFIX}.fdr0.01.hot.bed ${HD}/regions/${BASE}.broadPeak

done

