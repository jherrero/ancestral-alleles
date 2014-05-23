#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 SCRIPT fix_1kg_vcf.pl

This script has been used in May 2014 to fix the Ancestral Alleles information on the VCF files for the 1000 Genomes project.

Apart from using new AA calls for 1-bp indels, it also reformat the information such that the AA are now shown as an additional INFO tag,
suppressing the redundancy in the previous files.

The original 1000 Genomes VCFs files were compressed with tabix. The code uses gunzip to uncompress the files on the fly although it can work
with uncompressed files just as well.

The updated AA information was generated using the RunAncestralAllelesOnEMF.pm module. In order to reduce the number of predictions to be
calculated, that module can take a VCF file to limit the predictions to those positions. The aim of that exercise was to fix an issue with
the AA for 1-bp indels only. To achieve that efficiently, the 1-bp indels were extracted from the original VCF files and copied to a new set
of VCF files (one per chromosome).

=head2 1. Extracting 1-bp indels from the VCF files:

  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do \
    gunzip -c ../original_vcfs/ALL.chr$chr.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz | \
    perl -lane 'next if (/^#/); print if (($F[3] =~ /^[A-Z]$/ and $F[4] =~ /^[A-Z]{2}$/) or ($F[3] =~ /^[A-Z]{2}$/ and $F[4] =~ /^[A-Z]$/))' \
    > AA.chr$chr.indels.vcf; done

=head2 2. The files with the new predictions were then built using:

  standaloneJob.pl RunAncestralAllelesOnEMF.pm --emf emf/Compara.6_primates_EPO.chr1_1.emf.gz -vcf indel_vcfs/AA.chr1.indels.vcf -out AA.chr1_1.txt
  standaloneJob.pl RunAncestralAllelesOnEMF.pm --emf emf/Compara.6_primates_EPO.chr1_2.emf.gz -vcf indel_vcfs/AA.chr1.indels.vcf -out AA.chr1_2.txt
  standaloneJob.pl RunAncestralAllelesOnEMF.pm --emf emf/Compara.6_primates_EPO.chr1_3.emf.gz -vcf indel_vcfs/AA.chr1.indels.vcf -out AA.chr1_3.txt
  etc.

=head2 3. Finally, the new files were edited with this script:

  perl src/fix_1kg_vcf.pl original_vcfs/ALL.chr1.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz new_predictions/AA.chr1_* \
    | gzip -9 > edited_vcfs/ALL.chr1.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz
  perl src/fix_1kg_vcf.pl original_vcfs/ALL.chr2.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz new_predictions/AA.chr2_* \
    | gzip -9 > edited_vcfs/ALL.chr2.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz
  perl src/fix_1kg_vcf.pl original_vcfs/ALL.chr3.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz new_predictions/AA.chr3_* \
    | gzip -9 > edited_vcfs/ALL.chr3.vep_erb_gerp.20130502.all_var_types.sites.vcf.gz
  etc.

=head2 Software versions

=over

=item Ortheus

master branch downloaded on May 18 2014 from https://github.com/benedictpaten/ortheus

=item muscle

MUSCLE v3.8.31 by Robert C. Edgar

=item eHive

Commit 1054adef663b12d545ca3edd4ea35ef800648553 from https://github.com/jherrero/ensembl-hive, a fork from https://github.com/ensembl/ensembl-hive

=back

=head2 Acknowledgements

The authors acknowledge the use of the UCL Legion High Performance Computing facility, and associated services, in the completion of this work.

=cut


my $original_vcf_file = shift(@ARGV);
my @new_predictions_files = @ARGV;


our %event_type = (
1 => 'deletion',
2 => 'complex_deletion',
4 => 'insertion',
5 => 'complex_insertion',
7 => 'unsure');


my $new_predictions = {};
foreach my $this_new_predictions_file (@new_predictions_files) {
    read_new_predictions($this_new_predictions_file, $new_predictions);
#    print scalar(keys %{$new_predictions->{1}}), "\n";
}

if ($original_vcf_file =~ /\.gz/) {
    open(VCF, "gunzip -c $original_vcf_file |") or die;
} else {
    open(VCF, $original_vcf_file) or die;
}

while (<VCF>) {
    chomp;
    if ($_ =~ /^#/)
        {
        if ($_ eq '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|AA|IndelType|Ref|Alt">') {
            print '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF">',
                "\n",
            "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele, IndelType:Type of Indel (REF, ALT and IndelType are only defined to indels)\">\n";
        } else {
            print $_, "\n";
        }
        next;
        }
    
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split("\t", $_);

    my $edit_aa_info_for_this_vcf_entry = 0;
    my $event_type = undef;
    my $allele;
    if ($ref =~ /^[A-Z]$/ and $alt =~ /^[A-Z]{2}$/) {
        $edit_aa_info_for_this_vcf_entry = 1;
        $event_type = 'i';
        ($allele) = $alt =~ /^[A-Z]([A-Z])$/;
    }
    if ($ref =~ /^[A-Z]{2}$/ and $alt =~ /^[A-Z]$/) {
        $edit_aa_info_for_this_vcf_entry = 1;
        $event_type = 'd';
        ($allele) = $ref =~ /^[A-Z]([A-Z])$/;
    }
    
    my ($vep_info, $gerp1_info, $gerp2_info) = split(";", $info);
    my @consequences = split(",", $vep_info);
    my $new_info = "";
    my ($new_aa, $new_indel_type, $new_ref, $new_alt) = ("", "", "", "");
    my $new_ancestral_alleles_data_have_been_parsed = 0;
    foreach my $this_consequence (@consequences) {
        ## VCF_original: AA|IndelType|Ref|Alt
        ## $new_predictions: [IndelType, Ref, Alt, AA] or [Exception]
        ## VCF_new_format: Ref|Alt|AA|IndelType
        if (!$new_ancestral_alleles_data_have_been_parsed) {
            if ($edit_aa_info_for_this_vcf_entry) {
                if ($new_predictions->{$chr}->{$pos+1}->{$event_type}->{$allele}) {
                    ($new_indel_type, $new_ref, $new_alt, $new_aa) = @{$new_predictions->{$chr}->{$pos+1}->{$event_type}->{$allele}};
                    $new_indel_type = ($event_type{$new_indel_type} or $new_indel_type);
                } elsif ($new_predictions->{$chr}->{$pos+1}->{'e'}) {
                    my $exception = $new_predictions->{$chr}->{$pos+1}->{'e'};
                    if ($exception =~ /LOW_COMPLEXITY \(HR\)/) {
                        $new_indel_type = "unknown(HR)";
                    } elsif ($exception =~ /LOW_COMPLEXITY \(muscle (\d)/) {
                        $new_indel_type = "unknown(STR$1?)";
                    } else {
                        $new_indel_type = "unknown($exception)";
                    }
                } else {
                    $new_indel_type = "unknown(NO_COVERAGE)";
                }
            } else {
                ($new_aa) = $this_consequence =~ /\|([^\|]*)\|[^\|]*\|[^\|]*\|[^\|]*$/;
            }
            $new_ancestral_alleles_data_have_been_parsed = 1;
        }

        ## Remove current AA info from $this_consequence
        $this_consequence =~ s/\|[^\|]*\|[^\|]*\|[^\|]*\|[^\|]*$//;
        
        ## Add the resulting string to $new_info (which will be the final INFO string in the output)
        $new_info .= "," if ($new_info);
        $new_info .= $this_consequence;
    }

    $new_info .= ";AA=$new_aa|$new_ref|$new_alt|$new_indel_type" if ($new_aa or $new_ref or $new_alt or $new_indel_type);
    $new_info .= ";$gerp1_info" if ($gerp1_info);
    $new_info .= ";$gerp2_info" if ($gerp2_info);
#    print "$_\n";
    print join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $new_info), "\n";
    
}
close(VCF);

sub read_new_predictions {
    my ($file, $predictions) = @_;
    
    open(FILE, $file) or "die";
    while (<FILE>) {
        chomp;
        my @records = split(";", $_);
        my ($chr, $pos, $bp, $aa, $s) = split("\t", shift(@records));
#        print join(" -- ", $chr, $pos, $bp, $aa), "\n";
        foreach my $this_record (@records) {
            my @data = split("\t", $this_record);
            if (!@data) {
                # Nothing, formatting error for exceptions
            } elsif ($data[0] eq "e") {
                $predictions->{$chr}->{$pos}->{"e"} = $data[1];
            } elsif ($data[1] eq "i") {
                $predictions->{$chr}->{$pos}->{"i"}->{$data[0]} = [@data[2..5]];
            } elsif ($data[1] eq "d") {
                $predictions->{$chr}->{$pos}->{"d"}->{$data[0]} = [@data[2..5]];
            } else {
                die "Error: $_";
            }
        }
    }
    close(FILE);
}