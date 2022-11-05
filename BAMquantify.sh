#!/usr/bin/env bash
###
#Written by Luis Abatti <abatti.luis@gmail.com>, based on Linh Huynh's original script

# Function: Count reads at specific loci from .bed file for each bam file in the input directory (the first parameter below). Counts are then summarized into final_quantification.csv.
# Requirements
#    1. samtools version 1.10 or later
# Parameters
#    1. Path to the directory of bam files
#    2. Path to BED file containing coordinates (chr:start-end name)
#    3. SE or PE mode (SE/PE)
# Use
#    2. Run the script. For example: ./samtools_quantify_BAM.sh /usr/home/path_to_bam_files /usr/home/path_to_bed_file.bed PE
#    3. The result is summarized in the file final_quantification.csv
###

if ! [ -x "$(command -v samtools)" ]; then
  echo 'Error: samtools is not installed.' >&2
  exit 1
fi

ORIGIN_BAM_PATH=$1

if [[ -z "${ORIGIN_BAM_PATH}" ]]; then
    echo "Must provide BAM path directory" 1>&2
    exit 1
elif [[ -f "${ORIGIN_BAM_PATH}" ]]; then
    echo "BAM path is a file, not a directory" 1>&2
    exit 1
elif [[ -d "${ORIGIN_BAM_PATH}" ]]; then
    echo "Quantifying BAM files from ${ORIGIN_BAM_PATH}"
fi

BED_FILE=$2

if [[ -z "${BED_FILE}" ]]; then
    echo "Must provide .bed file" 1>&2
    exit 1
elif [[ ! -f "${BED_FILE}" ]]; then
    echo ".bed file does not exist" 1>&2
    exit 1
elif [[ -d "${BED_FILE}" ]]; then
    echo ".bed file path is a directory, not a file" 1>&2
    exit 1
elif [[ -f "${BED_FILE}" ]]; then
    echo "Quantifying loci from ${BED_FILE}"
fi

MODE=${3:-PE}

if [ "${MODE}" == "PE" ]
then
  echo "Paired-end mode"
elif [ "${MODE}" == "SE" ]
then
  echo "Single-end mode"
else
  echo "Please define MODE (PE/SE)" 1>&2
  exit 1
fi


chr_array=( $(awk '{print $1}' "${BED_FILE}" | sed -e 's/\r//g') )
start_array=( $(awk '{print $2}' "${BED_FILE}" | sed -e 's/\r//g') )
end_array=( $(awk '{print $3}' "${BED_FILE}" | sed -e 's/\r//g') )
name_array=( $(awk '{print $4}' "${BED_FILE}" | sed -e 's/\r//g') )
#sed removes the carriage return character

if [[ -z "${name_array}" ]]; then
  echo "${BED_FILE} does not contain a name column!" 1>&2
  exit 1
fi

OUTPUT_BAM_PATH="${ORIGIN_BAM_PATH}/Output_bam_files"

if [[ ! -d "${OUTPUT_BAM_PATH}" ]]; then
  mkdir "${OUTPUT_BAM_PATH}"
fi

output_file_list=""

for index in ${!name_array[*]}
#return the list of all array indices (0, 1, 2, 3...) then uses it to grab the position of each coordinate
  do
    out_name="${name_array[$index]}"
    output_bam_dir="${OUTPUT_BAM_PATH}/${out_name}"
    if [[ ! -d "${output_bam_dir}" ]]; then
      mkdir "${output_bam_dir}"
      #Makes directory if it doesn't exist
    fi
    pos="${index}"
    chr_pos="${chr_array[$pos]}"
    start_pos="${start_array[$pos]}"
    end_pos="${end_array[$pos]}"
    name_pos="${name_array[$pos]}"
    echo "Processing ${name_pos} ${chr_pos}:${start_pos}-${end_pos}"
    result_filename="${ORIGIN_BAM_PATH}/${out_name}.dat"
    if [[ -f ${result_filename} ]]; then
      rm "${result_filename}"
      #Removes result_filename if it already exists
    fi
    output_file_list="${output_file_list} ${result_filename}"
    if [ "${index}" -le 0 ]; then
      echo -e "file_name,read_count,${name_array[$index]}" >> "${result_filename}"
    else
      echo "${name_array[$index]}" >> "${result_filename}"
    fi
    for bam_file in $(ls "${ORIGIN_BAM_PATH}"/*.bam)
    do
      [[ -e "${bam_file}" ]] || break
      echo "Processing ${bam_file}"
      # remove the path
      tmp="${bam_file##*/}"
      # remove .bam
      filename="${tmp%%.bam}"
      bai_filename="${ORIGIN_BAM_PATH}/${filename}.bam.bai"
      output_bam_filename="${output_bam_dir}/${filename}.bam"
      if test -f "$bai_filename"; then
        echo "Index found: ${bai_filename}"
      else
        echo "Generating index: ${bai_filename}" 
        samtools index "${bam_file}" "${bai_filename}"
      fi
      if [ "${MODE}" == "PE" ]
        then
          samtools view -X -h "${bam_file}" "${bai_filename}" "${chr_pos}:${start_pos}-${end_pos}" > "${output_bam_filename}"
          read_num=$(samtools view -c -f 2 -q 10 "${output_bam_filename}")
            if [ "${index}" -le 0 ]; then
              total_read_num=$(samtools view -c -f 2 -q 10 "${bam_file}")
              echo -e "${filename},${total_read_num},${read_num}" >> "${result_filename}"
            else
              echo -e "${read_num}" >> "${result_filename}"
            fi
      elif [ "${MODE}" == "SE" ]
        then
          samtools view -X -h "${bam_file}" "${bai_filename}" "${chr_pos}:${start_pos}-${end_pos}" > "${output_bam_filename}"
          read_num=$(samtools view -c -q 10 "${output_bam_filename}")
            if [ "${index}" -le 0 ]; then
              total_read_num=$(samtools view -c -q 10 "${bam_file}")
              echo -e "${filename},${total_read_num},${read_num}" >> "${result_filename}"
            else
              echo -e "${read_num}" >> "${result_filename}"
            fi       
      fi 
    done
  done
paste -d "," ${output_file_list} > "${ORIGIN_BAM_PATH}"/final_quantification.csv
echo "Finished quantification!"
echo "Removing temporary files..."
rm ${output_file_list};
rm -r ${OUTPUT_BAM_PATH};
echo "Done!"
