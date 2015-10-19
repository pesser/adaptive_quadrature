if [ $# -ne 1 ]
then
  echo "Useage: $0 logfile"
  exit 1
fi
logfile=$1
grep -Eo "[0-9]*\.[0-9]* [0-9]*\.[0-9]* .*" $logfile > data
python visualize.py
