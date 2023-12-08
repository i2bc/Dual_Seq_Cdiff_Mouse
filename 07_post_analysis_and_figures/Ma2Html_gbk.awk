#!/bin/awk
# change gene annotation format
# from M2Html toward Genbank
# restrict to Strain CD630
# eg. :
# MaHtml geneID : CD1006 in column 2
# Genbank geneID : CD630_01006 
# Extract also functionnal categories (in column 1) and protein function (in column 2)
# eg. column 2:
# Sec-Translocase_Ribosomal_CD0058A-CD0060 
# => CD630_00580A, CD630_00590, CD630_00600
# ____________ BEGIN __________________
BEGIN{
   OFS="\t"
   FS="\t"
}
# __________ FUNCTIONS ________________
function letter_in_ID(gene_id, id){
# gene_id example: 0058A 0060
# id: 1=letter 2=number
#print "gene_id:",gene_id,"*"
   pos=match(gene_id,/[[:alpha:]]/)
   if(pos>0){
      id[1]=substr(gene_id,pos,length(gene_id))
      id[2]=substr(gene_id,1,pos-1)
      #print "letterinID:",id[1],"*gene_id:"gene_id"*"
   } else {
      id[1]=""
      id[2]=gene_id
   }
}
function print_with_zero(class,fun,gnumber,gletter){
   num=gnumber*1+0 # interger cast
   if (num<10){
      zero="000"
   } else { 
      if (num<100) {
         zero="00"
      } else {
         if (num<1000) {
            zero="0"
         } else {
            zero=""
         }
      }
   }
   if(fun~/_CD630/){
      if(num<10000){
         zero=zero"0";
      }
      printf("%s\t%s\tCD630_%s%d%s\n",class,fun,zero,gnumber,gletter)
   } else {   
      printf("%s\t%s\tCD630_%s%d0%s\n",class,fun,zero,gnumber,gletter)
   }
}
# _____________ MAIN ___________________
{
   if(($3=="CD630")&&($2!~/P0/)&&($2!~/P10/)){ # P0x or P10 = ID for plasmid sequence
      gsub(" ","_",$1)
      nbcol2=split($2,annot,"_")
      funct=annot[1]
      for (a=2;a<nbcol2;a++){
          funct=sprintf("%s_%s",funct,annot[a])
      } 
      gsub("CD","",annot[nbcol2])
      op=index(annot[nbcol2],"-") # last "-" == mark for operon
      if (op>0) { # operon case
         split(annot[nbcol2],op_ids,"-")
         letter_in_ID(op_ids[1], id_beg)
         letter_in_ID(op_ids[2], id_end)
         op_count=id_end[2]-id_beg[2] 
      } else { # ID with a unique gene
         letter_in_ID(annot[nbcol2], id_beg)
         op_count=0
      }
      print_with_zero($1,funct,id_beg[2]+0,id_beg[1]) # +0 for integer
      if (op>0) { # operon case
          for (i=1 ; i < op_count ; i++) {
             print_with_zero($1,funct,id_beg[2]+i,"")
          }
          print_with_zero($1,funct,id_end[2],id_end[1])
      }
   }
} 
