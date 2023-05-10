#cafetutorial_prep_r8s.py was downloaded from cafe turorial website
python cafetutorial_prep_r8s.py -i SpeciesTree_rooted.txt -o r8s_species11_file.txt -s 1000000 -p 'elegans_longest_transcript,briggsae_longest_transcript' -c '18.6'
r8s -b -f r8s_species11_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > twelve_spp_r8s_ultrametric.txt
awk -F"\t" '{print "Desc""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"' Orthogroups.GeneCount.tsv > cafe-input.data
sed -e '1 s/(null)/Desc/' -e '1 s/Orthogroup/Family/' -i cafe-input.data
python cafetutorial_clade_and_size_filter.py -i cafe-input.data -o filtered_cafe_input.txt -s 2> filtered.log
#cafe
mkdir -p reports
vi species11_cafe.sh
load -i filtered_cafe_input.txt -t 4 -l reports/log.txt -p 0.05
tree (((panamensis:11.951445,becei:11.951445)0.910146:10.310134,((inopinata:15.448491,elegans:15.448491)0.326161:3.151509,((((nigoni:4.208472,briggsae:4.208472)0.852374:7.602527,tribulationis:11.811000)0.540949:2.962477,(remanei:4.695755,latens:4.695755)0.839071:10.077722)0.352374:2.074682,tropicalis:16.848159)0.364371:1.751841)iptipt:3.661579)1:19.750842,bovis:42.012421)
lambda -s -t (((1,1)1,((1,1)1,((((1,1)1,1)1,(1,1)1)1,1)1)1)1,1)
report reports/report_run1_filter

chmod a+x species11_cafe.sh
cafe species11_cafe.sh
python cafetutorial_report_analysis.py -i reports/report_run1_filter.cafe -o reports/summary_run1
python cafetutorial_draw_tree.py -i reports/summary_run1_node.txt -t '(((panamensis:38.553206,becei:38.553206):33.258526,((inopinata:49.833869,elegans:49.833869):10.166131,((((nigoni:13.575658,briggsae:13.575658):24.524199,tribulationis:38.099857):9.556339,(remanei:15.147554,latens:15.147554):32.508642):6.692614,tropicalis:54.348810):5.651190):11.811732):63.713311,bovis:135.525043)' \
-d '(((panamensis<0>,becei<2>)<1>,((inopinata<4>,elegans<6>)<5>,((((nigoni<8>,briggsae<10>)<9>,tribulationis<12>)<11>,(remanei<14>,latens<16>)<15>)<13>,tropicalis<18>)<17>)<7>)<3>,bovis<20>)<19>' -o reports/summary_run1_tree_rapid.png -y Rapid

