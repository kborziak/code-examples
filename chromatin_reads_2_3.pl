#!/usr/bin/perl

use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# cut size = histone plus linker, ~200bp. sequence length is 51bp. +- 75 would center both plus and minus strand reads
#$fragment = 75;
$half = 100;
$fragment = 200;
$length = 51;
$chip_start = -2000;
$chip_end = 2000;
$window = 200;
$increment = 10;

$gff_file = "/home/kborziak/Documents/c_elegans/seeds/c_elegans.WS235.annotations.gff3";
$gene_names_file = "/home/kborziak/Documents/c_elegans/seeds/c_elegans.PRJNA13758.WS249.geneOtherIDs.txt";
$bin_file = '/home/kborziak/Documents/c_elegans/data/read_bins';
$bin_file2 = '/home/kborziak/Documents/c_elegans/data/reads_per_bin';

$files = {
#  'H3K4P' => '/home/kborziak/Documents/c_elegans/data/bam_parse_1_H3K4P', 
#  'H3K4C' => '/home/kborziak/Documents/c_elegans/data/bam_parse_1_H3K4C', 
#  'H3P' => '/home/kborziak/Documents/c_elegans/data/bam_parse_1_H3P', 
#  'H3C' => '/home/kborziak/Documents/c_elegans/data/bam_parse_1_H3C', 
  'wt_stv_con_1_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3K4C_1', 
  'wt_stv_con_1_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3C_1', 
  'wt_stv_con_2_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3K4C_2', 
  'wt_stv_con_2_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3C_2', 
  'wt_stv_pd_2_H3K4'   => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3K4P_1', 
  'wt_stv_pd_2_H3'     =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3P_1', 
  'wt_stv_pd_1_H3K4'   => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3K4P_2', 
  'wt_stv_pd_1_H3'     =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_NSH3P_2', 
  'csr_stv_pd_1_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPSH3K4_1', 
  'csr_stv_pd_1_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPSH3_1', 
  'csr_stv_pd_2_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPSH3K4_2', 
  'csr_stv_pd_2_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPSH3_2', 
  'csr_stv_con_2_H3K4' => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCSH3K4_1', 
  'csr_stv_con_2_H3'   =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCSH3_1', 
  'csr_stv_con_1_H3K4' => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCSH3K4_2', 
  'csr_stv_con_1_H3'   =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCSH3_2', 
  'csr_phe_pd_1_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPEH3K4_1', 
  'csr_phe_pd_1_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPEH3_1', 
  'csr_phe_pd_2_H3K4'  => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPEH3K4_2', 
  'csr_phe_pd_2_H3'    =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WPEH3_2', 
  'csr_phe_con_2_H3K4' => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCEH3K4_1', 
  'csr_phe_con_2_H3'   =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCEH3_1', 
  'csr_phe_con_1_H3K4' => '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCEH3K4_2', 
  'csr_phe_con_1_H3'   =>   '/home/kborziak/Documents/c_elegans/data/bam_data/bam_parse_2_WCEH3_2', 
};
$ratio = {
  'H3K4P' => 1, 
  'H3K4C' => 1, 
};
$control = {
  'csr_phe_con_1_H3K4' => 'csr_phe_con_1_H3', 
  'csr_phe_con_2_H3K4' => 'csr_phe_con_2_H3', 
  'csr_phe_pd_1_H3K4'  => 'csr_phe_pd_1_H3', 
  'csr_phe_pd_2_H3K4'  => 'csr_phe_pd_2_H3', 
  'csr_stv_con_1_H3K4' => 'csr_stv_con_1_H3', 
  'csr_stv_con_2_H3K4' => 'csr_stv_con_2_H3', 
  'csr_stv_pd_1_H3K4'  => 'csr_stv_pd_1_H3', 
  'csr_stv_pd_2_H3K4'  => 'csr_stv_pd_2_H3', 
  'wt_stv_con_1_H3K4'  => 'wt_stv_con_1_H3', 
  'wt_stv_con_2_H3K4'  => 'wt_stv_con_2_H3', 
  'wt_stv_pd_1_H3K4'   => 'wt_stv_pd_1_H3', 
  'wt_stv_pd_2_H3K4'   => 'wt_stv_pd_2_H3', 
};
$read_count = {
  'H3K4P' =>      2871625, 
  'H3P' =>       16430387, 
  'H3K4C' =>      4930338, 
  'H3C' =>       19769198, 
  'wt_stv_con_1_H3K4'  => 24035347, 
  'wt_stv_con_1_H3'    => 40887921, 
  'wt_stv_con_2_H3K4'  => 11671301, 
  'wt_stv_con_2_H3'    => 19481330, 
  'wt_stv_pd_2_H3K4'   => 25409597, 
  'wt_stv_pd_2_H3'     => 29020505, 
  'wt_stv_pd_1_H3K4'   => 48932444, 
  'wt_stv_pd_1_H3'     => 27735348, 
  'csr_stv_pd_1_H3K4'  => 3607812, 
  'csr_stv_pd_1_H3'    => 5364202, 
  'csr_stv_pd_2_H3K4'  => 29695921, 
  'csr_stv_pd_2_H3'    => 29258742, 
  'csr_stv_con_2_H3K4' => 37881131, 
  'csr_stv_con_2_H3'   => 42752980, 
  'csr_stv_con_1_H3K4' =>  2675327, 
  'csr_stv_con_1_H3'   =>  5879745, 
  'csr_phe_pd_1_H3K4'  => 22117611, 
  'csr_phe_pd_1_H3'    => 42566157, 
  'csr_phe_pd_2_H3K4'  => 47224258, 
  'csr_phe_pd_2_H3'    => 42976484, 
  'csr_phe_con_2_H3K4' => 22341783, 
  'csr_phe_con_2_H3'   => 42723514, 
  'csr_phe_con_1_H3K4' => 38375031, 
  'csr_phe_con_1_H3'   => 60148456, 
};
$chr_size = {
  'I' => 15072434, 
  'II' => 15279421, 
  'III' => 13783801, 
  'IV' => 17493829, 
  'V' => 20924180, 
  'X' => 17718942, 
};

open out_bins, ">", $bin_file;
open out_k4, ">", $bin_file."_k4_norm";
open out_h3, ">", $bin_file."_h3_norm";

open fin_norm, ">", $bin_file2;
open fin_k4, ">", $bin_file2."_k4_norm";


foreach $key (keys $control) {
  $coverage = 0; $overlap = 0;
  $out_file = '/home/kborziak/Documents/c_elegans/data/bam_data/'.$key.'_bed';
  open out, ">", $out_file;
  print out "track name=$key useScore=1\n";
  
  foreach $chr (keys $chr_size) {
    print "running $chr\n";
    undef $data;
    
    print "reading $key\n";
    open in, "<", $files->{$key};
    while ($line = <in>) {
      chomp $line;
      @split = split (/\t/, $line);
      if ($split[0] eq $chr) {
        if ($split[2] ne '') {
          $bin = sprintf ("%d", $split[1]/1000);
          $data->{$bin}->{$key} += $split[2];
        }
        if ($split[3] ne '') {
          $bin = sprintf ("%d", ($split[1] - $fragment + $length)/1000);
          $data->{$bin}->{$key} += $split[2];
        }
      }
    }
    close in;
    
    print "reading ".$control->{$key}."\n";
    open in, "<", $files->{$control->{$key}};
    while ($line = <in>) {
      chomp $line;
      @split = split (/\t/, $line);
      if ($split[0] eq $chr) {
        if ($split[2] ne '') {
          $bin = sprintf ("%d", $split[1]/1000);
          $data->{$bin}->{$control->{$key}} += $split[2];
        }
        if ($split[3] ne '') {
          $bin = sprintf ("%d", ($split[1] - $fragment + $length)/1000);
          $data->{$bin}->{$control->{$key}} += $split[2];
        }
      }
    }
    close in;
  
#  foreach $chr (keys $chr_size) {
    print "generating windows\n";
    $start = '';
    for ($i = 0; $i <= $chr_size->{$chr} /1000; $i++) {
      if (defined $data->{$i}->{$key} && $data->{$i}->{$key}  != 0) {
        $val2 = max ($data->{$i}->{$key} / $read_count->{$key} * 10**8 * (1 - ($data->{$i}->{$control->{$key}} / $read_count->{$control->{$key}} * 10**8) / ($data->{$i}->{$key} / $read_count->{$key} * 10**8)), 0);
      }
      else {
        $val2 = 0;
      }
      $round = sprintf ("%.2f", $val2);
      $k4 = sprintf ("%.2f", $data->{$i}->{$key} / $read_count->{$key} * 10**8);
      $h3 = sprintf ("%.2f", $data->{$i}->{$control->{$key}} / $read_count->{$control->{$key}} * 10**8);
#      $k4 = sprintf ("%.2f", $data->{$i}->{$key});
#      $h3 = sprintf ("%.2f", $data->{$i}->{$control->{$key}});
      $out_data->{'norm'}->{$round}->{$key}++;
      $out_data->{'k4'}->{$k4}->{$key}++;
      $out_data->{'h3'}->{$h3}->{$key}++;
      $fin_data->{'norm'}->{$chr}->{$i}->{$key} = $round;
      $fin_data->{'k4'}->{$chr}->{$i}->{$key} = $k4;
      $start = $i *1000;
      $end = min ($chr_size->{$chr}, $start + 999);
      print out "chr$chr\t$start\t$end\tval\t$val2\n";
      $start = '';
    }
  }
  close out;
}

print out_bins "bin";
foreach $key (keys $control) {
  print out_bins "\t$key";
}
print out_bins "\n";
foreach $val (sort {$a <=> $b} keys $out_data->{'norm'}) {
  print out_bins $val;
  foreach $key (keys $control) {
    if (defined $out_data->{'norm'}->{$val}->{$key}) {
      print out_bins "\t".$out_data->{'norm'}->{$val}->{$key};
    }
    else {
      print out_bins "\t0";
    }
  }
  print out_bins "\n";
}
close out_bins;

print out_k4 "bin";
foreach $key (keys $control) {
  print out_k4 "\t$key";
}
print out_k4 "\n";
foreach $val (sort {$a <=> $b} keys $out_data->{'k4'}) {
  print out_k4 $val;
  foreach $key (keys $control) {
    if (defined $out_data->{'k4'}->{$val}->{$key}) {
      print out_k4 "\t".$out_data->{'k4'}->{$val}->{$key};
    }
    else {
      print out_k4 "\t0";
    }
  }
  print out_k4 "\n";
}
close out_k4;

print out_h3 "bin";
foreach $key (keys $control) {
  print out_h3 "\t$key";
}
print out_h3 "\n";
foreach $val (sort {$a <=> $b} keys $out_data->{'h3'}) {
  print out_h3 $val;
  foreach $key (keys $control) {
    if (defined $out_data->{'h3'}->{$val}->{$key}) {
      print out_h3 "\t".$out_data->{'h3'}->{$val}->{$key};
    }
    else {
      print out_h3 "\t0";
    }
  }
  print out_h3 "\n";
}
close out_h3;


print fin_norm "chr\tbin";
foreach $key (keys $control) {
  print fin_norm "\t$key";
}
print fin_norm "\n";
foreach $chr (sort {$a <=> $b} keys $fin_data->{'norm'}) {
  foreach $bin (sort {$a <=> $b} keys $fin_data->{'norm'}->{$chr}) {
    print fin_norm "$chr\t$bin";
    foreach $key (keys $control) {
      if (defined $fin_data->{'norm'}->{$chr}->{$bin}->{$key}) {
        print fin_norm "\t".$fin_data->{'norm'}->{$chr}->{$bin}->{$key};
      }
      else {
        print fin_norm "\t0";
      }
    }
    print fin_norm "\n";
  }
}
close fin_norm;

print fin_k4 "chr\tbin";
foreach $key (keys $control) {
  print fin_k4 "\t$key";
}
print fin_k4 "\n";
foreach $chr (sort {$a <=> $b} keys $fin_data->{'k4'}) {
  foreach $bin (sort {$a <=> $b} keys $fin_data->{'k4'}->{$chr}) {
    print fin_k4 "$chr\t$bin";
    foreach $key (keys $control) {
      if (defined $fin_data->{'k4'}->{$chr}->{$bin}->{$key}) {
        print fin_k4 "\t".$fin_data->{'k4'}->{$chr}->{$bin}->{$key};
      }
      else {
        print fin_k4 "\t0";
      }
    }
    print fin_k4 "\n";
  }
}
close fin_k4;

