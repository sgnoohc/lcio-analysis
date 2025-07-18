if [ -z "$1" ]; then
      echo "Error: No argument supplied. Usage: $0 <RUNTAG>"
        exit 1
fi

RUNTAG="$1"

mkdir -p data/${RUNTAG}/

cp -r data/template/* data/${RUNTAG}/

OUTPUT_FILE="data/${RUNTAG}/run_sim_git_diffs.txt"
> "$OUTPUT_FILE"  # Truncate or create the output file

for dir in ../*/ ; do
  if [ -d "$dir/.git" ]; then
    echo "=== $dir ===" >> "$OUTPUT_FILE"
    (cd "$dir" && git diff) >> "$OUTPUT_FILE"
    echo -e "\n" >> "$OUTPUT_FILE"
  fi
done

ddsim \
    --steeringFile ../SteeringMacros/Sim/sim_steer_GEN_CONDOR.py \
    --inputFiles data/${RUNTAG}/muonGun_pT_0_50_gen.slcio \
    --numberOfEvents 1000 \
    --skipNEvents 0 \
    --outputFile data/${RUNTAG}/muonGun_pT_0_50_sim_0.slcio \
    --compactFile ../detector-simulation/geometries/MAIA_v0/MAIA_v0.xml
