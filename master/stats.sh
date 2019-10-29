a="Mesh_3d.java test_Mesh_3d.m  Tree_3d.java display_3d.m Tree_1d.java  Key2.java Tree.java Tree_2d.java test_dist3d.java Dist.java Key3.java test_object.java pdeeig_3d.m Potential_3D_Coulomb.class assemble_3d_phi_phi_java.m Potential_3D.java assemble_3d_dphi_dphi_java.m"


echo $a | wc -w 

# lines total
cat $a | wc -l 

# without blank lines
cat $a | grep -v "^\s*$"  | wc -l

# with neither blank lines nor comment lines
cat $a | grep -v "^.*//.*$" |grep -v "^\s*$"  | wc -l

#grep "(" tags $a | wc -l

#grep $a | grep -v "^.*//.*$" |grep -v "^\s*$"  | wc -l
