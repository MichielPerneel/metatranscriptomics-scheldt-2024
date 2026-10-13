#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
SRC="/user/scratch/gent/vo/000/gvo00081/vsc47816/Metatranscriptome_Scheldt_2025_05"
DEST="/user/scratch/gent/vo/000/gvo00081/vsc47816/metatranscriptomics-workflow/data/raw"

mkdir -p "$DEST"

# Collect: key = sample|R1 or sample|R2 ; value = list of files
declare -A files
 
while IFS= read -r -d '' f; do
  base="$(basename "$f")"              # N33_MKRN..._L8_1.fq.gz
  sample="${base%%_*}"                 # N33
  if [[ "$base" =~ _1\.fq\.gz$ ]]; then
    readtag="R1"
  elif [[ "$base" =~ _2\.fq\.gz$ ]]; then
    readtag="R2"
  else
    continue
  fi
 
  key="${sample}|${readtag}"
  files["$key"]+=$'\n'"$f"

done < <(find "$SRC" -type f -name "*.fq.gz" -print0)
 
# Process each sample/read group
for key in "${!files[@]}"; do
  sample="${key%%|*}"
  readtag="${key##*|}"
  out="$DEST/${sample}_${readtag}.fastq.gz"
 
  # Turn list into array, sort (keeps lane order stable)
  mapfile -t arr < <(printf "%s\n" "${files[$key]}" | sed '/^$/d' | sort)
 
  if (( ${#arr[@]} == 1 )); then
    # Single lane: symlink
    ln -sf "${arr[0]}" "$out"
    echo "LINK  $out -> ${arr[0]}"

  else
    # Multi-lane: concatenate into a new file
    tmp="${out}.tmp.$$"
    echo "MERGE $out  (${#arr[@]} files)"
    cat "${arr[@]}" > "$tmp"
    mv -f "$tmp" "$out"
  fi
done

echo "Done. Outputs in: $DEST"