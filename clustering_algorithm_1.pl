#!/usr/bin/perl 

use Math::Round qw (:all);
use List::Member;
use List::Util qw(min max sum shuffle);
use Scalar::Util qw(looks_like_number);
#use Data::Dumper;
use feature 'say';

sub choose {
  my $count = 0;
  my $result = undef;
  for( @_ ) {
    $result = $_ if rand() < 1 / ++$count;
  }
  $result;
}

$simulations = 1;
$mask_retrogenes = 0;
$intervening_genes = 2;
$cluster_size = 3;
$fraction_upregulated = 0;
$cluster_density = .66;
$max_cluster_size = 0;  # 0 for whole chromosome
$min_cluster_genes = 3;
$exp_level = 1;
$shuffle_types = {'none' => '', };#'all' => '', 'rgenes_only' => '', 'all_but_rgenes' => ''};
#$shuffle_used = 'none';# all or rgenes_only or all_but_rgenes or none
$algorithm_types = {'GLCA' => '', 'DCCA' => ''};
$algorithm = 'GLCA';# GLCA or DCCA
$output_clusters = 'yes';# yes or no

$tag = "NE_up_all";

print "Reading input\n";
$time1 = time;

#Clustering input
open elaine_in, "<", "/home/kborziak/Documents/c_elegans/data/gene_order2";
while ($elaine_line = <elaine_in>) {
  chomp ($elaine_line);
  @elaine_split = split ("\t", $elaine_line);
  if (looks_like_number($elaine_split[0])) {
    $ids->{$elaine_split[0]} = $elaine_split[2];
#    if ($elaine_split[10] > 0 && $elaine_split[7] ne 1) {
    if ($elaine_split[8] > 0) {
#    if ($elaine_split[7] eq 1) {
      $flyatlas_e->{$elaine_split[2]} = 1;
    }
    else {
      $flyatlas_e->{$elaine_split[2]} = 0;
    }
    if (!defined $chromosome_min->{$elaine_split[1]}) {
      $chromosome_min->{$elaine_split[1]} = $elaine_split[0];
    }
    $chromosome_max->{$elaine_split[1]} = $elaine_split[0];
    $chromosome->{$elaine_split[0]} = $elaine_split[1];
  }
}
close elaine_in;

#Pearson's correlation data
open pearson_file, "<", "/home/kborziak/Documents/retrogene/Pearsons_cor_t";
while ($pearson_line = <pearson_file>) {
  chomp ($pearson_line);
  @pearson_split = split ("\t", $pearson_line);
  $pearson_dat->{"$pearson_split[0] $pearson_split[1]"} = 0;
  $pearson_dat->{"$pearson_split[1] $pearson_split[0]"} = 0;
}
close pearson_file;

$time2 = time - $time1;
print "elapsed $time2 seconds\n";

#foreach $algorithm (keys %$algorithm_types) {
  foreach $shuffle_used (keys %$shuffle_types) {

#DCCA sim cluster definition & shuffle
if ($algorithm eq 'DCCA') {
  print "Starting DCCA shuffle $shuffle_used sim\n";
  open comp_out, ">", "/home/kborziak/Documents/retrogene/DCCA_sim_shuffle_".$shuffle_used."_3";
  print comp_out "simulation_#\t#_testis_clusters\ttestis_cluster_focal_genes\ttestis_genes_in_testis_clusters\tgenes_in_testis_clusters\t#_other_clusters\tother_cluster_focal_genes\ttestis_genes_in_other_clusters\tgenes_in_other_clusters\ttestis_retrogenes_in_testis_clusters\tother_retrogenes_in_testis_clusters\ttestis_retrogenes_in_other_clusters\tother_retrogenes_in_other_clusters\tfocal_testis_retrogenes_in_testis_clusters\tfocal_other_retrogenes_in_testis_clusters\tfocal_testis_retrogenes_in_other_clusters\tfocal_other_retrogenes_in_other_clusters";
#  for ($i = $cluster_size; $i <= 50; $i++) {
#    for ($i2 = $cluster_size; $i2 <= $i; $i2++) {
#      print comp_out "\ttestis_".$i."_focal_".$i2;#\tall_$i";
#      $cluster_stat_w_total->{$i}->{$i2}->{"test"} = 0;
#    }
#  }
#  for ($i = $cluster_size; $i <= 50; $i++) {
#    for ($i2 = $cluster_size; $i2 <= $i; $i2++) {
#      print comp_out "\tother_".$i."_focal_".$i2;#\tall_$i";
#      $cluster_stat_w_total->{$i}->{$i2}->{"other"} = 0;
#    }
#  }
  print comp_out "\tlarger_clusters\n";
  $time1 = time;
  for ($num = 1; $num <= $simulations; $num++) {
    $j = 0; $total = keys %$flyatlas_id; @shuffle = 1..$total; @shuffle = shuffle @shuffle; $b = 1;
    undef $ids2; undef $flyatlas_id2; undef $used;
    undef $cluster_stat_w; undef $cluster_stat_w_2; undef $cluster_stat_w_focal; undef $rgene_stat_w; undef $tgene_stat_w; undef $cluster_testis; undef $cluster_all; undef $cluster_focal;

#Shuffle types
    foreach $i (@shuffle) {
      if ($shuffle_used eq 'rgenes_only') {
        if (defined ($rgenes->{$ids2->{$i}})) {
          $ids2->{$b} = $ids->{$i};
          $id_shuff->{$b} = $i;
          $flyatlas_id2->{$b} = $flyatlas_id->{$i};
          $ids2->{$i} = $ids->{$b};
          $id_shuff->{$i} = $b;
          $flyatlas_id2->{$i} = $flyatlas_id->{$b};
        }
        elsif (!defined ($ids2->{$b})) {
          $ids2->{$b} = $ids->{$b};
          $id_shuff->{$b} = $b;
          $flyatlas_id2->{$b} = $flyatlas_id->{$b};
        }
      }
      elsif ($shuffle_used eq 'all') {
        $ids2->{$b} = $ids->{$i};
        $id_shuff->{$b} = $i;
        $flyatlas_id2->{$b} = $flyatlas_id->{$i};
      }
      elsif ($shuffle_used eq 'all_but_rgenes') {
        if (defined ($rgenes->{$ids2->{$b}})) {
          $ids2->{$b} = $ids->{$b};
          $id_shuff->{$b} = $b;
          $flyatlas_id2->{$b} = $flyatlas_id->{$b};
          $ids3->{$i} = $ids->{$i};
          $id_shuff->{$i} = $i;
          $flyatlas_id2->{$i} = $flyatlas_id->{$i};
        }
        elsif (!defined ($ids2->{$b})) {
          $ids2->{$b} = $ids->{$i};
          $id_shuff->{$b} = $i;
          $flyatlas_id2->{$b} = $flyatlas_id->{$i};
        }
      }
      elsif ($shuffle_used eq 'none') {
        $ids2->{$i} = $ids->{$i};
        $id_shuff->{$i} = $i;
        $flyatlas_id2->{$i} = $flyatlas_id->{$i};
      }
      $b++;
    }

#Dynamic clustering
    for ($i = 1; $i < $total - 1; $i++) {
      $x = 'go'; $temp = $i; $k = 0; @cluster = (); $exp = 0;
      if ($mask_retrogenes == 1 && defined ($rgenes->{$ids2->{$temp}})) {
        $used->{$temp} = 0;
      }
      elsif (!defined ($used->{$temp}) && $ids2->{$temp} ne '') {
        $exp = 0;
        while ($x eq 'go') {
          $temp = $temp + $k;
          $used->{$temp} = 0;
          push (@cluster, $temp);
          $k = 0; $l = 0;
          if ($flyatlas_e->{$ids2->{$temp}} >= $exp_level) {
            $exp++;
          }
          while ($l <= $intervening_genes) {
            $n = 0; $tot1 = 0; $tot2 = 0; $k++; $l++;
            if ($mask_retrogenes == 1 && defined ($rgenes->{$ids2->{$temp + $k}})) {
              $l--;
            }
            else {
              if (defined ($flyatlas_id2->{$temp}) && defined ($flyatlas_id2->{$temp + $k}) && !defined ($used->{$temp + $k})) {
                foreach $probe_1 (keys $flyatlas_id2->{$temp}) {
                  foreach $probe_2 (keys $flyatlas_id2->{$temp + $k}) {
                    if (defined ($pearson_dat->{"$probe_1 $probe_2"})) {
                      $n++;
                    }
                  }
                }
                $tot1 = keys $flyatlas_id2->{$temp};
                $tot2 = keys $flyatlas_id2->{$temp + $k};
              }
              if ($n == $tot1 * $tot2 && $n != 0) {
                $x = 'go';
                last;
              }
              $x = 'stop';
            }
          }
        }
      }
      if (@cluster >= $cluster_size) {
        $j++; $size = @cluster; $size2 = max (@cluster) - min (@cluster) + 1; $density = $exp / $size; $density2 = $exp / $size2; $tr_sim = 0; $or_sim = 0; $tr_sim_f = 0; $or_sim_f = 0; $t_gene = 0;
        foreach $id (@cluster) {
          $cluster_focal->{$id} = $cluster_focal->{$id}."cluster_$j,";
        }
        for ($id = min (@cluster); $id <= max (@cluster); $id++) {
          $cluster_all->{$id} = $cluster_all->{$id}."cluster_$j,";
          if ($flyatlas_e->{$ids2->{$id}} >= $exp_level) {
            $t_gene++;
          }
          if (defined ($rgenes->{$ids2->{$id}})) {
			      if ($id ~~ @cluster) {
              if ($flyatlas_e->{$ids2->{$id}} >= $exp_level) {
                $tr_sim_f++;
              }
              else {
                $or_sim_f++;
              }
            }
            else {
              if ($flyatlas_e->{$ids2->{$id}} >= $exp_level) {
                $tr_sim++;
              }
              else {
                $or_sim++;
              }
            }
          }
        }
        if ($density >= .75) {
          $clust_type = 'test';
          for ($id = min (@cluster); $id <= max (@cluster); $id++) {
            $cluster_testis->{$id} = $cluster_testis->{$id}."cluster_$j,";
          }
        }
        else {
          $clust_type = 'other';
        }
        $cluster_stat_w->{$size2}->{$size}->{$clust_type}++;
        $cluster_stat_w_2->{$size2}->{$clust_type}++;
        $cluster_stat_w_focal->{$size}->{$clust_type}++;
        $tgene_stat_w->{$size2}->{$clust_type} = $tgene_stat_w->{$size2}->{$clust_type} + $t_gene;
        $rgene_stat_w->{$size2}->{$clust_type."_test"} = $rgene_stat_w->{$size2}->{$clust_type."_test"} + $tr_sim;
        $rgene_stat_w->{$size2}->{$clust_type."_other"} = $rgene_stat_w->{$size2}->{$clust_type."_other"} + $or_sim;
        $rgene_stat_w->{$size2}->{$clust_type."_test_f"} = $rgene_stat_w->{$size2}->{$clust_type."_test_f"} + $tr_sim_f;
        $rgene_stat_w->{$size2}->{$clust_type."_other_f"} = $rgene_stat_w->{$size2}->{$clust_type."_other_f"} + $or_sim_f;
      }
    }
    
#Simulation output
    $tc_tg = 0; $oc_tg = 0; $tc_ag = 0; $oc_ag = 0; $t_c = 0; $o_c = 0; $tc_tr = 0; $oc_tr = 0; $tc_or = 0; $oc_or = 0; $tc_tr_f = 0; $oc_tr_f = 0; $tc_or_f = 0; $oc_or_f = 0; $tc_f = 0; $oc_f = 0;
    for ($key = $cluster_size; $key <= 50; $key++) {
      if (defined $cluster_stat_w->{$key}) {
        $tc_tg = $tc_tg + $tgene_stat_w->{$key}->{"test"};
        $oc_tg = $oc_tg + $tgene_stat_w->{$key}->{"other"};
        $tc_ag = $tc_ag + $key * $cluster_stat_w_2->{$key}->{"test"};
        $oc_ag = $oc_ag + $key * $cluster_stat_w_2->{$key}->{"other"};
        $t_c = $t_c + $cluster_stat_w_2->{$key}->{"test"};
        $o_c = $o_c + $cluster_stat_w_2->{$key}->{"other"};
        $tc_f = $tc_f + $key * $cluster_stat_w_focal->{$key}->{"test"};
        $oc_f = $oc_f + $key * $cluster_stat_w_focal->{$key}->{"other"};;
        $tc_tr = $tc_tr + $rgene_stat_w->{$key}->{"test_test"};
        $tc_or = $tc_or + $rgene_stat_w->{$key}->{"test_other"};
        $oc_tr = $oc_tr + $rgene_stat_w->{$key}->{"other_test"};
        $oc_or = $oc_or + $rgene_stat_w->{$key}->{"other_other"};
        $tc_tr_f = $tc_tr_f + $rgene_stat_w->{$key}->{"test_test_f"};
        $tc_or_f = $tc_or_f + $rgene_stat_w->{$key}->{"test_other_f"};
        $oc_tr_f = $oc_tr_f + $rgene_stat_w->{$key}->{"other_test_f"};
        $oc_or_f = $oc_or_f + $rgene_stat_w->{$key}->{"other_other_f"};
      }
    }
    print comp_out "$num\t$t_c\t$tc_f\t$tc_tg\t$tc_ag\t$o_c\t$oc_f\t$oc_tg\t$oc_ag\t"."$tc_tr\t$tc_or\t$oc_tr\t$oc_or\t$tc_tr_f\t$tc_or_f\t$oc_tr_f\t$oc_or_f";
    for ($key = $cluster_size; $key <= 50; $key++) {
      for ($i2 = $cluster_size; $i2 <= $key; $i2++) {
        if (defined $cluster_stat_w->{$key}->{$i2}->{"test"}) {
#          print comp_out "\t".$cluster_stat_w->{$key}->{$i2}->{"test"};#."\t".0;#.$cluster_stat_e->{$key}->{"other"};
          $cluster_stat_w_total->{$key}->{$i2}->{"test"} += $cluster_stat_w->{$key}->{$i2}->{"test"};
        }
        else {
#          print comp_out "\t0";#\t0";
        }
      }
    }
    for ($key = $cluster_size; $key <= 50; $key++) {
      for ($i2 = $cluster_size; $i2 <= $key; $i2++) {
        if (defined $cluster_stat_w->{$key}->{$i2}->{"other"}) {
#          print comp_out "\t".$cluster_stat_w->{$key}->{$i2}->{"other"};
          $cluster_stat_w_total->{$key}->{$i2}->{"other"} += $cluster_stat_w->{$key}->{$i2}->{"other"};
        }
        else {
#          print comp_out "\t0";
        }
      }
    }
    print comp_out "\t";
    foreach $key (keys %$cluster_stat_w) {
      if ($key > 50) {
        for ($i2 = $cluster_size; $i2 <= $key; $i2++) {
          if (defined $cluster_stat_w->{$key}->{$i2}) {
#            print comp_out $key.": ".$i2.": ";
            if (defined $cluster_stat_w->{$key}->{$i2}->{"test"}) {
#              print comp_out $cluster_stat_w->{$key}->{$i2}->{"test"};
              $cluster_stat_w_total->{$key}->{$i2}->{"test"} += $cluster_stat_w->{$key}->{$i2}->{"test"};
            }
#            print comp_out ",";
            if (defined $cluster_stat_w->{$key}->{$i2}->{"other"}) {
#              print comp_out $cluster_stat_w->{$key}->{$i2}->{"other"};
              $cluster_stat_w_total->{$key}->{$i2}->{"other"} += $cluster_stat_w->{$key}->{$i2}->{"other"};
            }
#            print comp_out "; ";
          }
        }
      }
    }
    print comp_out "\n";
  }
  close comp_out;
  $time2 = time - $time1;
  print "elapsed $time2 seconds\n";
  if ($output_clusters eq 'yes') {
    open out_clusters, ">", "/home/kborziak/Documents/retrogene/".$algorithm."_sim_suffle_".$shuffle_used."_clusters_2";
    print out_clusters "id\tgene_name\tfocal_genes\ttestis_clusters\tall_clusters\n";
    for ($i = 1; $i <= $total; $i++) {
      print out_clusters $i."\t".$ids->{$i}."\t".$cluster_focal->{$i}."\t".$cluster_testis->{$i}."\t".$cluster_all->{$i}."\n";
    }
#    print out_clusters "size\tfocal\ttestis_count\tother_count\n";
    foreach $key (sort { $a <=> $b} keys %$cluster_stat_w_total) {
      foreach $i2 (sort { $a <=> $b} keys $cluster_stat_w_total->{$key}) {
#        print out_clusters $key."\t".$i2."\t";
        if (defined $cluster_stat_w_total->{$key}->{$i2}->{"test"}) {
#          print out_clusters $cluster_stat_w_total->{$key}->{$i2}->{"test"};
        }
#        print out_clusters "\t";
        if (defined $cluster_stat_w_total->{$key}->{$i2}->{"other"}) {
          print out_clusters $cluster_stat_w_total->{$key}->{$i2}->{"other"};
        }
#        print out_clusters "\n";
      }
    }
    close out_clusters;
  }
}

#GLCA sim cluster definition & shuffle
elsif ($algorithm eq 'GLCA') {
  print "Starting GLCA shuffle $shuffle_used sim\n";
  open comp_out, ">", "/home/kborziak/Documents/c_elegans/data/$algorithm"."_".$tag."_sim_shuffle_".$shuffle_used."_2";
  print comp_out "simulation_#\t#_clusters\t#_focal_clusters\ttotal_genes_in_clusters";
  for ($i = $cluster_size; $i <= 50; $i++) {
    for ($i2 = .5; $i2 <= 1; $i2 += .1) {
      print comp_out "\ttestis_".$i."_focal_".$i2;#\tall_$i";
      $cluster_stat_e_total->{$i}->{$i2}->{"test"} = 0;
    }
  }
  for ($i = $cluster_size; $i <= 50; $i++) {
    for ($i2 = $cluster_size; $i2 <= $i; $i2++) {
#      print comp_out "\tother_".$i."_focal_".$i2;#\tall_$i";
      $cluster_stat_e_total->{$i}->{$i2}->{"other"} = 0;
    }
  }
  print comp_out "\tlarger_clusters\n";
  $time1 = time;
  for ($num = 1; $num <= $simulations; $num++) {
    $j = 1; $total = keys %$ids; @shuffle = 1..$total; @shuffle = shuffle @shuffle; $b = 1;
    undef $ids3; undef $flyatlas_id3; undef $cluster_esim;
    undef $cluster_stat_e; undef $cluster_stat_e_2; undef $rgene_stat_e; undef $tgene_stat_e; undef $cluster_testis; undef $cluster_all; undef $cluster_focal;

#Shuffle types
    foreach $i (@shuffle) {
      if ($shuffle_used eq 'rgenes_only') {
        if (defined ($rgenes->{$ids3->{$i}})) {
          $ids3->{$b} = $ids->{$i};
          $id_shuff->{$b} = $i;
          $flyatlas_id3->{$b} = $flyatlas_id->{$i};
          $ids3->{$i} = $ids->{$b};
          $id_shuff->{$i} = $b;
          $flyatlas_id3->{$i} = $flyatlas_id->{$b};
        }
        elsif (!defined ($ids3->{$b})) {
          $ids3->{$b} = $ids->{$b};
          $id_shuff->{$b} = $b;
          $flyatlas_id3->{$b} = $flyatlas_id->{$b};
        }
      }
      elsif ($shuffle_used eq 'all') {
        $ids3->{$b} = $ids->{$i};
        $id_shuff->{$b} = $i;
#        $flyatlas_id3->{$b} = $flyatlas_id->{$i};
      }
      elsif ($shuffle_used eq 'all_but_rgenes') {
        if (defined ($rgenes->{$ids3->{$b}})) {
          $ids3->{$b} = $ids->{$b};
          $id_shuff->{$b} = $b;
          $flyatlas_id3->{$b} = $flyatlas_id->{$b};
          $ids3->{$i} = $ids->{$i};
          $id_shuff->{$i} = $i;
          $flyatlas_id3->{$i} = $flyatlas_id->{$i};
        }
        elsif (!defined ($ids3->{$b})) {
          $ids3->{$b} = $ids->{$i};
          $id_shuff->{$b} = $i;
          $flyatlas_id3->{$b} = $flyatlas_id->{$i};
        }
      }
      elsif ($shuffle_used eq 'none') {
        $ids3->{$i} = $ids->{$i};
        $id_shuff->{$i} = $i;
#        $flyatlas_id3->{$i} = $flyatlas_id->{$i};
      }
      $b++;
    }

#Forward
    undef $used;
    for ($i = 1; $i < $total; $i++) {
      $exp = 0;
      if ($mask_retrogenes == 1 && !defined ($rgenes->{$ids3->{$i}})) {
        $used->{$i} = 0;
      }
      elsif (!defined ($used->{$i}) && $ids3->{$i} ne '') {
        if ($flyatlas_e->{$ids->{$id_shuff->{$i}}} >= $exp_level) {
          if ($max_cluster_size == 0) {
            $temp = $chromosome_max->{$chromosome->{$i}};
          }
          else {
            $temp = min ($i + $max_cluster_size, $chromosome_max->{$chromosome->{$i}});
          }
          while ($temp > $i) {
            if ($flyatlas_e->{$ids->{$id_shuff->{$temp}}} >= $exp_level) {
              $c_tot = $temp - $i + 1;
              for ($k = $i; $k <= $temp; $k++) {
                $exp = $exp + $flyatlas_e->{$ids3->{$k}};
              }
              if ($exp / $c_tot >= $cluster_density && $exp >= $min_cluster_genes) {
                for ($k = $i; $k <= $temp; $k++) {
                  $used->{$k} = 0;
                  if ($flyatlas_e->{$ids3->{$k}} == 1) {
                    $cluster_esim->{$k} = 1;
                  }
                  else {
                    $cluster_esim->{$k} = 0;
                  }
                }
                $temp = $i;
              }
              else {
                $exp = 0;
                $temp--;
              }
            }
            else {
              $temp--;
            }
          }
        }
      }
    }

#Reverse
    undef $used;
    for ($i = $total; $i > 1; $i--) {
      $exp = 0;
      if ($mask_retrogenes == 1 && !defined ($rgenes->{$ids3->{$temp}})) {
        $used->{$i} = 0;
      }
      elsif (!defined ($used->{$i}) && $ids3->{$i} ne '') {
        if ($flyatlas_e->{$ids->{$id_shuff->{$i}}} >= $exp_level) {
          if ($max_cluster_size == 0) {
            $temp = $chromosome_min->{$chromosome->{$i}};
          }
          else {
            $temp = max ($i - $max_cluster_size, $chromosome_min->{$chromosome->{$i}});
          }
          while ($temp < $i) {
            if ($flyatlas_e->{$ids->{$id_shuff->{$temp}}} >= $exp_level) {
              $c_tot = $i - $temp + 1;
              for ($k = $temp; $k <= $i; $k++) {
                $exp = $exp + $flyatlas_e->{$ids3->{$k}};
              }
              if ($exp / $c_tot >= $cluster_density && $exp >= $min_cluster_genes) {
                for ($k = $temp; $k <= $i; $k++) {
                  $used->{$k} = 0;
                  if ($flyatlas_e->{$ids3->{$k}} == 1) {
                    $cluster_esim->{$k} = 1;
                  }
                  else {
                    $cluster_esim->{$k} = 0;
                  }
                }
                $temp = $i;
              }
              else {
                $exp = 0;
                $temp++;
              }
            }
            else {
              $temp++;
            }
          }
        }
      }
    }
    
#Cluster definition
    if ($num eq 1 && $shuffle_used eq "none") {
      open comp_out_2, ">", "/home/kborziak/Documents/c_elegans/data/$algorithm"."_".$tag."_clusters";
      print comp_out_2 "Order\tChromosome\tID\tCluster\n";
    }
    $size = 0; $exp = 0; $tr_sim = 0;  $or_sim = 0; $tg_sim = 0; $focal_size = 0; $tr_sim_f = 0;  $or_sim_f = 0; $cluster_number = 1;
    for ($i = 1; $i <= $total; $i++) {
      if ($num eq 1 && $shuffle_used eq "none") {
        print comp_out_2 "$i\t".$chromosome->{$i}."\t".$ids->{$i}."\t";
      }
      if (defined $cluster_esim->{$i}) {
        $size++;
        $exp = $exp + $cluster_esim->{$i};
        if (defined ($rgenes->{$ids3->{$i}})) {
          if ($flyatlas_e->{$ids3->{$i}} >= $exp_level) {
            $tr_sim_f++;
          }
          else {
            $or_sim++;
          }
        }
        $cluster_testis->{$i} = $cluster_testis->{$i}."cluster_$j,";
        $cluster_all->{$i} = $cluster_all->{$i}."cluster_$j,";
        if ($flyatlas_e->{$ids3->{$i}} >= $exp_level) {
          $cluster_focal->{$i} = $cluster_focal->{$i}."cluster_$j,";
          $focal_size++;
        }
        if ($num eq 1 && $shuffle_used eq "none") {
          print comp_out_2 "cluster_$cluster_number";
        }
      }
      elsif ($size > 0) {
        $j++; $size_arr_e->{$size} = '';
        $round = $focal_size / $size;
        $round = sprintf "%.1f", $round;
        $cluster_stat_e->{$size}->{$round}->{"test"}++;
#        $cluster_stat_e->{$size}->{$focal_size}->{"test"}++;
        $cluster_stat_e_2->{$size}->{"test"}++;
        $tgene_stat_e->{$size}->{"test"} = $tgene_stat_e->{$size}->{"test"} + $exp;
        $rgene_stat_e->{$size}->{"test_test"} = $rgene_stat_e->{$size}->{"test_test"} + $tr_sim;
        $rgene_stat_e->{$size}->{"test_other"} = $rgene_stat_e->{$size}->{"test_other"} + $or_sim;
        $rgene_stat_e->{$size}->{"test_test_f"} = $rgene_stat_e->{$size}->{"test_test_f"} + $tr_sim_f;
        $rgene_stat_e->{$size}->{"test_other_f"} = $rgene_stat_e->{$size}->{"test_other_f"} + $or_sim_f;
        $size = 0; $exp = 0; $tr_sim = 0;  $or_sim = 0; $tr_sim_f = 0;  $or_sim_f = 0;  $focal_size = 0;
        $cluster_number++;
      }
      if ($num eq 1 && $shuffle_used eq "none") {
        print comp_out_2 "\n";
      }
    }
    
#Simulation output
    $tc_tg = 0; $oc_tg = 0; $tc_ag = 0; $oc_ag = 0; $t_c = 0; $o_c = 0; $tc_tr = 0; $oc_tr = 0; $tc_or = 0; $oc_or = 0; $tc_tr_f = 0; $oc_tr_f = 0; $tc_or_f = 0; $oc_or_f = 0; $tc_f = 0; $oc_f = 0;
    for ($key = $min_cluster_genes; $key <= max(keys %$cluster_stat_e); $key++) {
      if (defined $cluster_stat_e->{$key}) {
        $tc_tg = $tc_tg +        $tgene_stat_e->{$key}->{"test"};
        $tc_ag = $tc_ag + $key * $cluster_stat_e_2->{$key}->{"test"};
        $t_c = $t_c +        $cluster_stat_e_2->{$key}->{"test"};
        $tc_f = $tc_f + $tgene_stat_e->{$key}->{"test"};
        $tc_tr = $tc_tr +    $rgene_stat_e->{$key}->{"test_test"};
        $tc_or = $tc_or +    $rgene_stat_e->{$key}->{"test_other"};
        $tc_tr_f = $tc_tr_f + $rgene_stat_e->{$key}->{"test_test_f"};
        $tc_or_f = $tc_or_f + $rgene_stat_e->{$key}->{"test_other_f"};
      }
    }
    print comp_out "$num\t$t_c\t$tc_f\t$tc_ag";
    for ($key = $min_cluster_genes; $key <= 50; $key++) {
      for ($i2 = .5; $i2 <= 1; $i2 += .1) {
        if (defined $cluster_stat_e->{$key}->{$i2}->{"test"}) {
          print comp_out "\t".$cluster_stat_e->{$key}->{$i2}->{"test"};#."\t".0;#.$cluster_stat_e->{$key}->{"other"};
          $cluster_stat_e_total->{$key}->{$i2}->{"test"} += $cluster_stat_e->{$key}->{$i2}->{"test"};
        }
        else {
          print comp_out "\t0";#\t0";
        }
      }
    }
    for ($key = $min_cluster_genes; $key <= 50; $key++) {
      for ($i2 = $cluster_size; $i2 <= $key; $i2++) {
        if (defined $cluster_stat_e->{$key}->{$i2}->{"other"}) {
#          print comp_out "\t".$cluster_stat_e->{$key}->{$i2}->{"other"};
          $cluster_stat_e_total->{$key}->{$i2}->{"other"} += $cluster_stat_e->{$key}->{$i2}->{"other"};
        }
        else {
#          print comp_out "\t0";
        }
      }
    }
    print comp_out "\t";
    foreach $key (keys %$cluster_stat_e) {
      if ($key > 50) {
        for ($i2 = .5; $i2 <= 1; $i2 += .1) {
#        for ($i2 = $cluster_size; $i2 <= $key; $i2++) {
          if (defined $cluster_stat_e->{$key}->{$i2}) {
            print comp_out $key.": ".$i2.": ";
            if (defined $cluster_stat_e->{$key}->{$i2}->{"test"}) {
              print comp_out $cluster_stat_e->{$key}->{$i2}->{"test"};
              $cluster_stat_e_total->{$key}->{$i2}->{"test"} += $cluster_stat_e->{$key}->{$i2}->{"test"};
            }
            print comp_out ",";
            if (defined $cluster_stat_e->{$key}->{$i2}->{"other"}) {
              print comp_out $cluster_stat_e->{$key}->{$i2}->{"other"};
              $cluster_stat_e_total->{$key}->{$i2}->{"other"} += $cluster_stat_e->{$key}->{$i2}->{"other"};
            }
            print comp_out "; ";
          }
        }
      }
    }
    print comp_out "\n";
  }
  close comp_out;
  $time2 = time - $time1;
  print "elapsed $time2 seconds\n";
  print Dumper $cluster_stat_e_total;
  
  if ($output_clusters eq 'yes') {
#    open out_clusters, ">", "/home/kborziak/Documents/retrogene/".$algorithm."_sim_suffle_".$shuffle_used."_clusters_2";
#    print out_clusters "id\tgene_name\tfocal_genes\ttestis_clusters\tall_clusters\n";
    for ($i = 1; $i <= $total; $i++) {
#      print out_clusters $i."\t".$ids->{$i}."\t".$cluster_focal->{$i}."\t".$cluster_testis->{$i}."\t".$cluster_all->{$i}."\n";
    }
#    print out_clusters "size\tfocal\ttestis_count\tother_count\n";
    foreach $key (sort { $a <=> $b} keys %$cluster_stat_e_total) {
      foreach $i2 (sort { $a <=> $b} keys $cluster_stat_e_total->{$key}) {
#        print out_clusters $key."\t".$i2."\t";
        if (defined $cluster_stat_e_total->{$key}->{$i2}->{"test"}) {
#          print out_clusters $cluster_stat_e_total->{$key}->{$i2}->{"test"};
        }
#        print out_clusters "\t";
        if (defined $cluster_stat_e_total->{$key}->{$i2}->{"other"}) {
#          print out_clusters $cluster_stat_e_total->{$key}->{$i2}->{"other"};
        }
#        print out_clusters "\n";
      }
    }
#    close out_clusters;
  }
}

  }
#}

