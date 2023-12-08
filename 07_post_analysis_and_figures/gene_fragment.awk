#/bin/awk
# fragmented gene management 
# x=part number
# duplicate x times the ligne 
# add "_px" at the geneID end
#
BEGIN{
   nblist=split(listp,pp,",")
}
{
   if(index(listp,$3)>0){ # gene in fragment
      gene_notfind=1 
      i=1 
      while ((gene_notfind==1)&&(i<=nblist)){
         if(pp[i]~$3){
           gene_notfind=0
           nbp=pp[i+1]
           for(p=1;p<=nbp;p++){
              printf("%s\t%s\t%s_p%d\n",$1,$2,$3,p)
           }
        }
        i=i+2
     }
  } else { # gene without fragment
     print $0
  }
}
