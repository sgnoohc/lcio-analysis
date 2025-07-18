if [ -z "$1" ]; then
      echo "Error: No argument supplied. Usage: $0 <RUNTAG>"
        exit 1
fi

RUNTAG="$1"

sh run_sim.sh ${RUNTAG}
sh run_reco.sh ${RUNTAG}
