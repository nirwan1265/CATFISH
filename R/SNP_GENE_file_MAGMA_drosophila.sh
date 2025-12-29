# ### First download the gff file from flybase
# # https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#GFF_files

# # Gene loc:
# awk -F'\t' '
# BEGIN{OFS="\t"}
# $1 ~ /^(2L|2R|3L|3R|4|X|Y)$/ && $3=="gene" {
#   id="";
#   n=split($9,a,";");
#   for(i=1;i<=n;i++){
#     if(a[i] ~ /^ID=FBgn/){
#       sub(/^ID=/,"",a[i]);
#       id=a[i];
#       break
#     }
#   }
#   if(id!="") print id, $1, $4, $5
# }' dmel*.gff > dmel.flybase.fbgn.genes.loc


# # Mapping file
# awk -F'\t' '
# BEGIN{OFS="\t"; print "FBgn","symbol","CG_alias","chr","start","end"}
# $1 ~ /^(2L|2R|3L|3R|4|X|Y)$/ && $3=="gene" {
#   id=""; nm=""; alias=""; cg="";
#   n=split($9,a,";");
#   for(i=1;i<=n;i++){
#     if(a[i] ~ /^ID=FBgn/){ tmp=a[i]; sub(/^ID=/,"",tmp); id=tmp; }
#     if(a[i] ~ /^Name=/){  tmp=a[i]; sub(/^Name=/,"",tmp); nm=tmp; }
#     if(a[i] ~ /^Alias=/){ tmp=a[i]; sub(/^Alias=/,"",tmp); alias=tmp; }
#   }
#   # pick the first CGxxxxx from Alias list if present
#   if(alias!=""){
#     m=split(alias,b,",");
#     for(j=1;j<=m;j++){
#       if(b[j] ~ /^CG[0-9]+$/){ cg=b[j]; break }
#     }
#   }
#   if(id!="") print id, nm, cg, $1, $4, $5
# }' dmel*.gff > dmel.flybase.gene_map.tsv



# # R code for gene l