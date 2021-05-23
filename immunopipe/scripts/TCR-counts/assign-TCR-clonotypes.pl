#! /usr/bin/env perl


use IO::File;
use Getopt::Std;
getopts("o:s:p:");
$outdir = $opt_o;
$output = $opt_s;
@prefixes = split ",",$opt_p;

# $vdj_dir = "../Data/GSE139555";
# $expr_dir = "../Data/GSE139555";

# @samples = @ARGV;
$len_samples = $#ARGV - 1;
@vdj_samples = @ARGV[0..$len_samples/2];
@rna_samples = @ARGV[$len_samples/2+1..$len_samples+1];

# Learn the universe of clonotypes and print master list
$OUT = new IO::File(">$outdir/$output.master");
$clonotype_idx = 1;

foreach $sample (@vdj_samples) {
    $last_barcode = "";
    @alpha_nts = ();
    @beta_nts = ();
    %aa = ();

    $FP = new IO::File("cat $sample |") or die;
    if (defined($line = <$FP>)) {
	# Skip header
    }

    while (defined($line = <$FP>)) {
	chop $line;
	# print STDERR "Line is $line\n";
	@fields = split ",",$line;
	$barcode = $fields[0];
	$chain = $fields[5];
	$tenx_clonotype = $fields[16];
	$clonotype_aa = $fields[12];
	$clonotype_nt = $fields[13];

	if ($fields[17] !~ /consensus/) {
	    # Skip

	} elsif ($clonotype_nt eq "None") {
	    # Skip

	} elsif ($chain ne "TRA" && $chain ne "TRB") {
	    # Focus only on alpha and beta for now

	} elsif ($last_barcode eq "") {
	    # Initial barcode
	    $aa{$clonotype_nt} = $clonotype_aa;
	    if ($chain eq "TRA") {
		push @alpha_nts,$clonotype_nt;
	    } elsif ($chain eq "TRB") {
		push @beta_nts,$clonotype_nt;
	    }

	    $last_barcode = $barcode;
	    $last_tenx_clonotype = $tenx_clonotype;

	} elsif ($barcode eq $last_barcode) {
	    # Continuation of previous barcode
	    $aa{$clonotype_nt} = $clonotype_aa;
	    if ($chain eq "TRA") {
		push @alpha_nts,$clonotype_nt;
	    } elsif ($chain eq "TRB") {
		push @beta_nts,$clonotype_nt;
	    }

	} else {
	    # New barcode observed.  First, handle the previous barcode
	    # print STDERR "Barcode $barcode is different from last one $last_barcode\n";
	    @alpha_sorted = sort {$a cmp $b} @alpha_nts;
	    @beta_sorted = sort {$a cmp $b} @beta_nts;
	    $alpha_consensus = join(",",@alpha_sorted);
	    $beta_consensus = join(",",@beta_sorted);
	    if (defined($id = $clonotype_id{$alpha_consensus . "|" . $beta_consensus})) {
		# Existing universal clonotype.  Check for consistency with CellRanger clonotype assignments
		# print STDERR "Saw clonotype $id again with alpha consensus $alpha_consensus and beta consensus $beta_consensus\n";
		if (!defined($tenx_clonotype{$sample}{$id})) {
		    # Can happen with a clonotype from a previous sample
		    $tenx_clonotype{$sample}{$id} = $last_tenx_clonotype;
		} elsif ($tenx_clonotype{$sample}{$id} ne $last_tenx_clonotype) {
		    print STDERR "For sample $sample, clonotype $id, 10x clonotype has changed from $tenx_clonotype{$sample}{$id} to $last_tenx_clonotype\n";
		}

	    } else {
			# New universal clonotype
			$id = $output . ".C" . $clonotype_idx;
			$clonotype_id{$alpha_consensus . "|" . $beta_consensus} = $id;
			$tenx_clonotype{$sample}{$id} = $last_tenx_clonotype;
			# print STDERR "Assigning $id to $last_tenx_clonotype $alpha_consensus | $beta_consensus\n";
			print $OUT $id . "\t";
			print $OUT "$alpha_consensus\t";
			print $OUT join(",",map {$aa{$_}} @alpha_sorted) . "\t";
			print $OUT "$beta_consensus\t";
			print $OUT join(",",map {$aa{$_}} @beta_sorted) . "\n";
			push @clonotype_ids,$id;
			$clonotype_idx++;
	    }
	    $assignment{$sample}{$last_barcode} = $id;


	    # Now handle the first entry for this new barcode
	    @alpha_nts = ();
	    @beta_nts = ();
	    %aa = ();
	    $aa{$clonotype_nt} = $clonotype_aa;
	    if ($chain eq "TRA") {
		push @alpha_nts,$clonotype_nt;
	    } elsif ($chain eq "TRB") {
		push @beta_nts,$clonotype_nt;
	    }

	    $last_barcode = $barcode;
	    $last_tenx_clonotype = $tenx_clonotype;
	}
    }

    # print STDERR "End of file\n";
    @alpha_sorted = sort {$a cmp $b} @alpha_nts;
    @beta_sorted = sort {$a cmp $b} @beta_nts;
    $alpha_consensus = join(",",@alpha_sorted);
    $beta_consensus = join(",",@beta_sorted);
    if (defined($id = $clonotype_id{$alpha_consensus . "|" . $beta_consensus})) {
	# print STDERR "Saw clonotype $id again with alpha consensus $alpha_consensus and beta consensus $beta_consensus\n";
	if (!defined($tenx_clonotype{$sample}{$id})) {
	    # Can happen with a clonotype from a previous sample
	    $tenx_clonotype{$sample}{$id} = $last_tenx_clonotype;
	} elsif ($tenx_clonotype{$sample}{$id} ne $last_tenx_clonotype) {
	    print STDERR "For sample $sample, clonotype $id, 10x clonotype has changed from $tenx_clonotype{$sample}{$id} to $last_tenx_clonotype\n";
	}

    } else {
		$id = $output . ".C" . $clonotype_idx;
		$clonotype_id{$alpha_consensus . "|" . $beta_consensus} = $id;
		$tenx_clonotype{$sample}{$id} = $last_tenx_clonotype;
		# print STDERR "Assigning $id to $last_tenx_clonotype $alpha_consensus | $beta_consensus\n";
		print $OUT $id . "\t";
		print $OUT "$alpha_consensus\t";
		print $OUT join(",",map {$aa{$_}} @alpha_sorted) . "\t";
		print $OUT "$beta_consensus\t";
		print $OUT join(",",map {$aa{$_}} @beta_sorted) . "\n";
		push @clonotype_ids,$id;
		$clonotype_idx++;
    }
    $assignment{$sample}{$last_barcode} = $id;

    close($FP);
}
close($OUT);


@header = ("alpha.V","alpha.J","alpha.C","beta.V","beta.D","beta.J","beta.C","gamma","delta","igh","igk","igl","clonotype","clonotype.10x");
$nfields = $#header + 1;


# Generate phenotype table and tally counts
$samplei = 1;
foreach $sample (@vdj_samples) {
    $prefix = shift(@prefixes);

    # Gather
    $last_barcode = "";
    $last_tenx_clonotype = "";
    $tra_v_gene = $tra_d_gene = $tra_j_gene = $tra_c_gene = "NA";
    $trb_v_gene = $trb_d_gene = $trb_j_gene = $trb_c_gene = "NA";
    $trg_gene = $trd_gene = "NA";
    $igh_gene = $igk_gene = $igl_gene = "NA";
    $clonotype_id = "NA";

    $FP = new IO::File("cat $sample |") or die;
    if (defined($line = <$FP>)) {
	# Skip header
    }

    while (defined($line = <$FP>)) {
	chop $line;
	@fields = split ",",$line;
	$barcode = $fields[0];
	$tenx_clonotype = $fields[16];

	if ($barcode ne $last_barcode) {
	    if ($last_barcode ne "") {
		if (!defined($clonotype_id = $assignment{$sample}{$last_barcode})) {
		    # print STDERR "Sample $sample, barcode $last_barcode has no clonotype\n";
		    $clonotype_id = "NA";
		}
		if ($last_tenx_clonotype eq "") {
		    $last_tenx_clonotype = "NA";
		}
		@ {$clonotype_info{$sample}{$last_barcode}} = ($tra_v_gene,$tra_j_gene,$tra_c_gene,$trb_v_gene,$trb_d_gene,$trb_j_gene,$trb_c_gene,
							       $trg_gene,$trd_gene,
							       $igh_gene,$igk_gene,$igl_gene,$clonotype_id,$last_tenx_clonotype);
	    }
	    $tra_v_gene = $tra_d_gene = $tra_j_gene = $tra_c_gene = "NA";
	    $trb_v_gene = $trb_d_gene = $trb_j_gene = $trb_c_gene = "NA";
	    $trg_gene = $trd_gene = "NA";
	    $igh_gene = $igk_gene = $igl_gene = "NA";
	    $clonotype_id = "NA";
	    $last_barcode = $barcode;
	    $last_tenx_clonotype = $tenx_clonotype;
	}

	# Chains: TRA, TRB, TRG, TRD, IGH, IGK, IGL, Multi
	if (($chain = $fields[5]) eq "None") {
	    # Skip
	} elsif ($chain eq "TRA") {
	    if ($fields[6] ne "None") {
		$tra_v_gene = $fields[6];
	    }
	    if ($fields[7] ne "None") {
		$tra_d_gene = $fields[7];
	    }
	    if ($fields[8] ne "None") {
		$tra_j_gene = $fields[8];
	    }
	    if ($fields[9] ne "None") {
		$tra_c_gene = $fields[9];
	    }

	} elsif ($chain eq "TRB") {
	    if ($fields[6] ne "None") {
		$trb_v_gene = $fields[6];
	    }
	    if ($fields[7] ne "None") {
		$trb_d_gene = $fields[7];
	    }
	    if ($fields[8] ne "None") {
		$trb_j_gene = $fields[8];
	    }
	    if ($fields[9] ne "None") {
		$trb_c_gene = $fields[9];
	    }

	} elsif ($chain eq "TRG") {
	    $trg_gene = "TRG";
	    #if ($fields[8] ne "None") {
	    #	$trg_j_gene = $fields[8];
	    #}
	    #if ($fields[9] ne "None") {
	    #	$trg_c_gene = $fields[9];
	    #}

	} elsif ($chain eq "TRD") {
	    $trd_gene = "TRD";
	    #if ($fields[6] ne "None") {
	    #	$trd_v_gene = $fields[6];
	    #}
	    #if ($fields[8] ne "None") {
	    #	$trd_j_gene = $fields[6];
	    #}

	} elsif ($chain eq "IGH") {
	    $igh_gene = "IGH";

	} elsif ($chain eq "IGK") {
	    $igk_gene = "IGK";

	} elsif ($chain eq "IGL") {
	    $igl_gene = "IGL";

	} elsif ($chain eq "Multi") {
	    #if ($fields[6] =~ /TRD/) {
	    #	$trd_v_gene = $fields[6];
	    #}
	    #if ($fields[8] =~ /TRD/) {
	    #	$trd_j_gene = $fields[6];
	    #}
	}

    }
    close($FP);

    if ($last_barcode ne "") {
	if (!defined($clonotype_id = $assignment{$sample}{$last_barcode})) {
	    # print STDERR "Sample $sample, barcode $last_barcode has no clonotype\n";
	    $clonotype_id = "NA";
	} else {
	}
	if ($last_tenx_clonotype eq "") {
	    $last_tenx_clonotype = "NA";
	}
	@ {$clonotype_info{$sample}{$last_barcode}} = ($tra_v_gene,$tra_j_gene,$tra_c_gene,$trb_v_gene,$trb_d_gene,$trb_j_gene,$trb_c_gene,
						       $trg_gene,$trd_gene,$igh_gene,$igk_gene,$igl_gene,$clonotype_id,$last_tenx_clonotype);
    }


    # Print
    $OUT = new IO::File(">$outdir/$output.pheno$samplei");
    print $OUT "barcode" . "\t";
    print $OUT join("\t",@header) . "\n";

    # Only 1% of raw samples have a clonotype, but 80% of the filtered ones do

    # $file = "$expr_dir/$sample.barcodes.tsv.gz";
    $file = @rna_samples[$samplei-1];
    $FP = new IO::File("gunzip -f -c $file |") or warn "No barcodes.tsv file for $sample.  Cannot combine with expression data";
    if (!defined($FP)) {
	foreach $barcode (keys % {$clonotype_info{$sample}}) {
	    print $OUT $prefix . "_" . $barcode . "\t";   # Add prefix to match barcodes in Seurat
	    print $OUT join("\t",@ {$clonotype_info{$sample}{$barcode}}) . "\n";
	}

    } else {
	while (defined($line = <$FP>)) {
	    chop $line;
	    ($barcode) = $line =~ /(\S+)/;

	    print $OUT $prefix . "_" . $barcode . "\t";   # Add prefix to match barcodes in Seurat

	    if (!defined($clonotype_info{$sample}{$barcode})) {
		print $OUT "NA";
		for ($i = 1; $i < $nfields; $i++) {
		    print $OUT "\t" . "NA"
		}
		print $OUT "\n";

	    } else {
		print $OUT join("\t",@ {$clonotype_info{$sample}{$barcode}}) . "\n";
		$clonotype_id = $ {$clonotype_info{$sample}{$barcode}}[12];
		$ {$counts[$samplei]}{$clonotype_id} += 1;
	    }
	}
	close($FP);
    }

    close($OUT);

    $samplei++;
}

########################################################################
#   Counts
########################################################################

$OUT = new IO::File(">$outdir/$output.counts");
foreach $clonotype_id (@clonotype_ids) {
    print $OUT $clonotype_id;
    for ($samplei = 1; $samplei <= $#samples + 1; $samplei++) {
	if (!defined($count = $ {$counts[$samplei]}{$clonotype_id})) {
	    print $OUT "\t" . "0";
	} else {
	    print $OUT "\t" . $count;
	}
    }
    print $OUT "\n";
}
close($OUT);


exit;
