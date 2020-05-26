
while read ID
do
file=$(ls $SOURCE/*.vcf | grep $ID)
if [[ $file == "" ]];
then
continue
fi
grep "frameshift" $file | sed "s/frameshift_variant/ /g" | 
while read line; 
do 
protein=$(echo $line | awk '{print $4}' | awk -F'|' '{print $3}'); 
echo -e "$ID\t$protein" >> frameshift.genelist; 
done; 
done < "$INPUT"
