#!/usr/bin/env bash
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":m:e:p:s:h:" arg; do
  case $arg in
    m) # Specify count matrix.
      matrix=${OPTARG}
      echo "m is ${OPTARG}"
      ;;
    e) # Specify extended count matrix (incl gene lengths)
      extmatrix=${OPTARG}
      echo "e is ${OPTARG}"
      ;;
    p) # Specify prefix name for output files
      prefix=${OPTARG}
      echo "p is ${OPTARG}"
      ;;
    s) # Specify file with sample data (incl. sequencing depth)
      sampledata=${OPTARG}
      echo "s is ${OPTARG}"
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

#Run the scripts
Rscript /home/people/pmun/bin/MGAS/countmatrixinfo2FPKM.R -m ${extmatrix} -s ${sampledata} -p ${prefix}
Rscript /home/people/pmun/bin/MGAS/countmatrix2diversity.R -m ${matrix} -p ${prefix}
Rscript /home/people/pmun/bin/MGAS/countmatrix2metaMDS.R -m ${matrix} -p ${prefix}
Rscript /home/people/pmun/bin/MGAS/matrix2heatmap.R -m ${prefix}_FPKM_matrix.txt -p ${prefix}
