if test $# -lt 1; then
   echo Usage:  ./caver_rename.sh caver_folder
   exit 0
fi

for i in $1/*; do 
   j=$1/$(basename $i | sed 's/.*\.\([0-9]*\)$/\1/').pdb
   echo Changing file $i to $j... 
   mv $i $j
done

