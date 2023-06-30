if test $# -lt 1; then
   echo "Usage: ./check_log.sh logfile"
   exit 0
fi

if sed '/nohup/'d $1 | grep -q '^[^0-9]'; then
   echo 
   echo 'Warning: potential issues detected'
   echo 
   echo '===================='
   sed '/nohup/'d $1 | grep -A 1 -B 1 '^[^0-9]'
   echo 
   echo 
fi

 
