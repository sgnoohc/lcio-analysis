if [ -z "$1" ]; then
      echo "Error: No argument supplied. Usage: $0 <RUNTAG>"
        exit 1
fi

RUNTAG="$1"
export RUNTAG=${RUNTAG}

mkdir -p data/${RUNTAG}/

OUTPUT_FILE="data/${RUNTAG}/run_reco_git_diffs.txt"
> "$OUTPUT_FILE"  # Truncate or create the output file

for dir in ../*/ ; do
  if [ -d "$dir/.git" ]; then
    echo "=== $dir ===" >> "$OUTPUT_FILE"
    (cd "$dir" && git diff) >> "$OUTPUT_FILE"
    echo -e "\n" >> "$OUTPUT_FILE"
  fi
done

k4run \
    ../SteeringMacros/k4Reco/steer_reco.py \
    --TypeEvent muonGun_pT_0_50
