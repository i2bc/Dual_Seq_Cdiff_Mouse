#
# this awk script gets header line and adds 2 columns
# run with:
# awk -f thisAwkScript.awk inputFile.txt > outputFile.txt
#
BEGIN{
  FS="\t"
}
{
  if(NR==1){
    print(sprintf("%s\tcogClass\tprod",$0))
  }
}
