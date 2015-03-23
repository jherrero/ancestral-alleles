#! /usr/bin/env perl
use strict;
use warnings;

use BaseAncestralAlleles;
use Getopt::Long;

my $help;
my $emf;
my $species = "homo_sapiens";
my $maf;
my $input;
my $output;
my $format = "tsv";
my $ortheus_exe = "/home/regmher/src/OrtheusC_2010-01-18/bin/OrtheusC";
my $muscle_exe = "~/src/muscle3.8.31_i86linux64";
# my $ortheus_exe = "ortheus_core";
# my $muscle_exe = "muscle";
my $work_dir = "/tmp";
my $max_alignment_length = 100;
my $flank_length = 20;
    
my $verbose = 0;
my $alignments_file = "";
my $allow_gzip;

my $EXCEPTION;

my $desc = qq{
USE: run_aa.pl [options] --emf EMF_DIR --input INPUT.txt --output OUTPUT.txt

OPTIONS:
 --help
    Shows this help.

 --emf EMF_DIR
    The location of the EMF directory with an index file (homo_sapiens.index or similar for other
    species).

 --species SPECIES
    The name of the reference species.
    [Def: $species] (other species untested)
 
 --input INPUT.txt
    The name of the input file. See --format option
    Supports compressed files (.gz or .bz2)
 
 --format tsv
    The only format currently support is tsv.
    [Def: $format]
 
 --output OUTPUT.txt
    The name of the output file.
    Supports compressed files (.gz or .bz2)
    [Def: }.($output?$output:"STDOUT").qq{]
 
 --work-dir TEMP_DIR
    The name of the temporary directory
    [Def: $work_dir]
 
 --allow-gzip
    Allows to used compressed EMF files. THIS IS MUCH SLOWER THAN UNCOMPRESSING ALL THE EMF FILES
    PRIOR TO RUNNING THIS TOOL.

 --verbose
    Output the alignments to the standard output or to the ALIGNMENT.out file (see --alignment-file
    option).
    
 --alignment-file ALIGNMENT.out
    Stores the alignments in this file.
    Supports compressed files (.gz or .bz2)
    [Def: Don't store them unless verbose option is used]
 
 --ortheus ortheus
    Path for ortheus.
    [Def: $ortheus_exe]

 --muscle muscle
    Path for muscle.
    [Def: $muscle_exe]

 --flank-length $flank_length
    Sets the length of the flank in the reference sequence to extract and re-align the alleles.
    [Def: $flank_length]

 --max-alignment-length $max_alignment_length
    Refuse to infer the AA if the alignment for the flanking region is longer than this.
    [Def: $max_alignment_length]
};


GetOptions(
    "help" => \$help,
    "emf=s" => \$emf,
    "species=s" => \$species,
    "format=s" => \$format,
    "input=s" => \$input,
    "output=s" => \$output,
    "work-dir|work_dir=s" => \$work_dir,
    "allow-gzip|allow_gzip" => \$allow_gzip,
    "alignments-file|alignments_file=s" => \$alignments_file,
    "ortheus=s" => \$ortheus_exe,
    "muscle=s" => \$muscle_exe,
    "max-alignment-length|max_alignment_length=i" => \$max_alignment_length,
    "flank-length|flank_length=i" => \$flank_length,
    "verbose" => \$verbose,
    );

if (!$emf or !$input or $help) {
    print $desc;
    exit !$help;
}

my $emf_index_cache = get_emf_index_cache("$emf/$species.index");

if ($input =~ /\.bz2$/) {
    open(INPUT, "bunzip2 -c $input |") or die;
} elsif ($input =~ /\.gz$/) {
    open(INPUT, "gunzip -c $input |") or die;
} else {
    open(INPUT, $input) or die;
}

if ($output and $output =~ /\.bz2$/) {
    open(OUTPUT, "| bzip2 -c > $output") or die;
} elsif ($output and $output =~ /\.gz$/) {
    open(OUTPUT, "| gzip -c > $output") or die;
} elsif ($output) {
    open(OUTPUT, ">$output") or die;
} else {
    *OUTPUT = *STDOUT;
}

if ($alignments_file =~ /\.bz2$/) {
    $verbose = 1;
    open(VERBOSE, "| bzip2 -c > $alignments_file") or die;
} elsif ($alignments_file =~ /\.gz$/) {
    $verbose = 1;
    open(VERBOSE, "| gzip -c > $alignments_file") or die;
} elsif ($alignments_file) {
    $verbose = 1;
    open(VERBOSE, "> $alignments_file") or die;
} else {
    *VERBOSE = *STDOUT;
}

my $count = 0;
#$_ = <INPUT>;
#chomp;
#my @headers = split("\t", $_);
#print OUTPUT join("\t", @headers, "AA.call", "Alleles"), "\n";

while (<INPUT>) {
    my ($chr, $start, $end, $ref_allele, $alt_allele);
    chomp;
    my @data = split("\t", $_);
    if ($format eq "tsv") {
        $chr = $data[2];
        $start = $data[3];
        $end = $data[4];
        $ref_allele = $data[5];
        $alt_allele = $data[6];
    } else {
        die "Unknown format [$format]";
    }
    my $variant_type = ($ref_allele =~ tr/A-Za-z/A-Za-z/)."/".($alt_allele =~ tr/A-Za-z/A-Za-z/);

    print VERBOSE "INPUT: $_\n" if ($verbose);

    my ($ref, $alt, $anc, $call) = get_ancestral_allele($emf, $emf_index_cache, $chr, $start,
            $end, $ref_allele, $alt_allele);
    if (!$call) {
        if ($EXCEPTION) {
            print VERBOSE "EXCEPTION: $EXCEPTION\t$variant_type\n" if ($verbose);
            print OUTPUT join("\t", @data, $EXCEPTION, "$ref_allele/$alt_allele/?"), "\n";
            $EXCEPTION = undef;
        }
        next;
    }
    if ($verbose) {
        print VERBOSE "ANCESTRAL ALLELE: $ref/$alt/$anc ", join("/", map {length($_)} ($ref, $alt, $anc)), "\n";
        print VERBOSE "ANCESTRY CALL FOR ALT: $call\t$variant_type\t$ref/$alt/$anc ", join("/",
                map {length($_)} ($ref, $alt, $anc)), "\n\n";
    }

    print OUTPUT join("\t", @data, $call, "$ref/$alt/$anc"), "\n";
}

if ($output) {
    close(OUTPUT);
}
if ($alignments_file) {
    close(VERBOSE);
}

sub get_ancestral_allele {
    my ($emf, $emf_index_cache, $chr, $start, $end, $ref_allele, $alt_allele) = @_;
    my $ancestral_allele;
    (my $sub_sorted_alignment, $alt_allele) = get_emf_sub_alignment($emf, $emf_index_cache, $chr, $start, $end, $ref_allele, $alt_allele);
    return undef if (!$sub_sorted_alignment);
    
    ## Adjust for homopolymers
    

    my $reference_ortheus_alignment = run_ortheus($sub_sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose);
    get_alleles_from_ortheus_alignment($reference_ortheus_alignment, $flank_length, $flank_length);
    if ($verbose) {
        print VERBOSE "===============================================\n";
        print VERBOSE  " ALIGNMENT FOR REFERENCE ALLELE\n";
        print VERBOSE  "===============================================\n";
        print_ortheus_alignment($reference_ortheus_alignment, *VERBOSE);
        print VERBOSE  "===============================================\n";
        print VERBOSE  "\n";
    }

    $ref_allele =~ s/\-//g;
    $alt_allele =~ s/\-//g;
    my $alt_seq = $sub_sorted_alignment->{sequences}->[0]->{original_sequence};
#     print $alt_seq, "\n";
    $alt_seq =~ s/(.{$flank_length}).*(.{$flank_length})/$1$alt_allele$2/;
#     print $alt_seq, "\n";

    $sub_sorted_alignment->{sequences}->[0]->{original_sequence} = $alt_seq;
    my $alternate_ortheus_alignment = run_ortheus($sub_sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose);
    get_alleles_from_ortheus_alignment($alternate_ortheus_alignment, $flank_length, $flank_length);
    if ($verbose) {
        print VERBOSE  "===============================================\n";
        print VERBOSE  " ALIGNMENT FOR ALTERNATE ALLELE\n";
        print VERBOSE  "===============================================\n";
        print_ortheus_alignment($alternate_ortheus_alignment, *VERBOSE);
        print VERBOSE  "===============================================\n";
        print VERBOSE  "\n";
    }
    my ($ref, $alt, $anc, $call) = call_ancestral_allele_from_ortheus_alignments($reference_ortheus_alignment, $alternate_ortheus_alignment, $flank_length, $flank_length);

    return ($ref, $alt, $anc, $call);
}

sub get_emf_index_cache {
    my ($emf_index) = @_;
    my $emf_index_cache;

    die "Cannot find the EMF index file $emf_index\n" if (!-e $emf_index);
    open(INDEX, $emf_index) or die;
    my ($chr, $start, $end, $strand, $this_emf_file, $offset);
    while (<INDEX>) {
      chomp;
      ($chr, $start, $end, $strand, $this_emf_file, $offset) = split("\t", $_);
      push(@{$emf_index_cache->{$chr}}, [$start, $end, $strand, $this_emf_file, $offset]);
    }
    close(INDEX);

    foreach my $chr (keys %$emf_index_cache) {
        $emf_index_cache->{$chr} = [sort {$a->[0] <=> $b->[0]} @{$emf_index_cache->{$chr}}]
    }

    return $emf_index_cache;
}

sub get_emf_sub_alignment {
    my ($emf, $emf_index_cache, $chr, $start, $end, $ref_allele, $alt_allele) = @_;
    
    my $blocks = $emf_index_cache->{$chr};
    return undef if (!$blocks);
    my $min_i = 0;
    my $max_i = @$blocks - 1;
    my $i = int(($min_i + $max_i) / 2);
    my $index;
    for (; 1; ) {
        if ($max_i - $min_i < 5) {
            for ($i = $min_i; $i < $max_i; $i++) {
                if ($blocks->[$i]->[0] <= $start and $start <= $blocks->[$i]->[1]) {
                    $index = $i;
                    last;
                }
            }
            last;
        } elsif ($blocks->[$i]->[0] > $start) {
            $max_i = $i;
            $i = int(($min_i + $max_i) / 2);
        } elsif ($blocks->[$i]->[1] < $start) {
            $min_i = $i;
            $i = int(($min_i + $max_i) / 2);
        } else {
            $index = $i;
            last;
        }
    }
    if (!defined($index)) {
        $EXCEPTION = "NO.ALIGNMENT";
        return undef;
    } elsif ($blocks->[$i]->[0] > $start - $flank_length or $end + $flank_length > $blocks->[$i]->[1]) {
        $EXCEPTION = "NO.FLANK.ALIGNMENT";
        return undef;
#     } else {
#         print join(" - ", $blocks->[$i]->[0], $start - $flank_length, $end + $flank_length, $blocks->[$i]->[1]), "\n";
    }

    my $this_emf_file = $blocks->[$index]->[3];
    my $offset = $blocks->[$index]->[4];
    my $sorted_alignment = get_emf_original_alignment($this_emf_file, $offset);

    ## Reverse-complement if required
    straighten_sorted_alignment($sorted_alignment);
    
    my $ref_original_seq = $sorted_alignment->{sequences}->[0]->{aligned_sequence};
    $ref_original_seq =~ s/\-//g;
    ## Check ref allele
    if ($ref_allele eq "-" and $start == $end + 1) {
        # No sequence to check, ref allele is pad anyway.
    } elsif (uc(substr($ref_original_seq, $start - $blocks->[$i]->[0], $end - $start + 1)) ne uc($ref_allele)) {
        $EXCEPTION = "WRONG.REF.ALLELE (".uc(substr($ref_original_seq, $start - $blocks->[$i]->[0], $end - $start + 1)).")";
        return undef;
    }

    ## Re-define start and end for indels on homopolymer runs
    if ($ref_allele eq "-" and ($alt_allele =~ /^A+$/ or $alt_allele =~ /^C+$/ or $alt_allele =~ /^G+$/ or $alt_allele =~ /^T+$/)) {
        my $indel = substr($alt_allele, 0, 1);
        my ($extend5) = substr($ref_original_seq, 0, $start - $blocks->[$i]->[0]) =~ /($indel+)$/;
#         print "...", substr(substr($ref_original_seq, 0, $start - $blocks->[$i]->[0]), -10), "\n";
        if ($extend5) {
#             print "EXTEND5 $extend5\n";
            $start -= length($extend5);
            $alt_allele = $extend5.$alt_allele;
        }
        my ($extend3) = substr($ref_original_seq, $end - $blocks->[$i]->[0] + 1) =~ /^($indel+)/;
#         print substr(substr($ref_original_seq, $end - $blocks->[$i]->[0] + 1), 0, 10), "...\n";
        if ($extend3) {
            $end += length($extend3);
#             print "EXTEND3 $extend3\n";
            $alt_allele .= $extend3;
        }
        if ($blocks->[$i]->[0] > $start - $flank_length or $end + $flank_length > $blocks->[$i]->[1]) {
            $EXCEPTION = "NO.FLANK.ALIGNMENT";
            return undef;
        }
    }
    if ($alt_allele eq "-" and ($ref_allele =~ /^A+$/ or $ref_allele =~ /^C+$/ or $ref_allele =~ /^G+$/ or $ref_allele =~ /^T+$/)) {
        my $indel = substr($ref_allele, 0, 1);
        my ($extend5) = substr($ref_original_seq, 0, $start - $blocks->[$i]->[0]) =~ /($indel+)$/;
#         print "...", substr(substr($ref_original_seq, 0, $start - $blocks->[$i]->[0]), -10), "\n";
        if ($extend5) {
#             print "EXTEND5 $extend5\n";
            $start -= length($extend5);
            $alt_allele = $extend5.$alt_allele;
        }
        my ($extend3) = substr($ref_original_seq, $end - $blocks->[$i]->[0] + 1) =~ /^($indel+)/;
#         print substr(substr($ref_original_seq, $end - $blocks->[$i]->[0] + 1), 0, 10), "...\n";
        if ($extend3) {
            $end += length($extend3);
#             print "EXTEND3 $extend3\n";
            $alt_allele .= $extend3;
        }
        if ($blocks->[$i]->[0] > $start - $flank_length or $end + $flank_length > $blocks->[$i]->[1]) {
            $EXCEPTION = "NO.FLANK.ALIGNMENT";
            return undef;
        }
    }
        

    ## Calculate the length of the region to be trimmed from the 5' of the alignment and the length
    ## of the region to be kept ($trim_length_in_cols and length($aligned_allele))
    my $ref_aligned_seq = $sorted_alignment->{sequences}->[0]->{aligned_sequence};
    my $trim_length_in_cols = 0;
    my $trim_length_in_bp = $start - $sorted_alignment->{sequences}->[0]->{start} - $flank_length;
    my $allele_length_in_bp = $end - $start + 1;

    foreach my $x (100000, 10000, 1000) {
        while ($trim_length_in_bp > $x) {
            my $first_kb = substr($ref_aligned_seq, 0, $x, "");
            $trim_length_in_cols += $x;
            my $num_nucl = $first_kb =~ tr/ACTGNactgn/ACTGNactgn/;
            $trim_length_in_bp -= $num_nucl;
        }
    }
    my ($trim) = $ref_aligned_seq =~ /^(\-*(?:[a-zA-Z]\-*){$trim_length_in_bp})/;
    substr($ref_aligned_seq, 0, length($trim), "");
    $trim_length_in_cols += length($trim);
    my ($flank5) = $ref_aligned_seq =~ /^(\-*(?:[a-zA-Z]\-*){$flank_length})/;
    substr($ref_aligned_seq, 0, length($flank5), "");
    my ($aligned_allele) = $ref_aligned_seq =~ /^((?:\-*[a-zA-Z]){$allele_length_in_bp})/;
    substr($ref_aligned_seq, 0, length($aligned_allele), "");
    my ($flank3) = $ref_aligned_seq =~ /^((?:\-*[a-zA-Z]){$flank_length})/;

    my $alignment_length = length($flank5) + length($aligned_allele) + length($flank3);
#     print join(" :: ", $flank5, $aligned_allele, $flank3), "\n";


    my $reference_sequence_exception = get_reference_sequence_exception($flank5, $flank3, $alignment_length, $max_alignment_length, $muscle_exe, $work_dir, $verbose);
    if ($reference_sequence_exception) {
        $EXCEPTION = $reference_sequence_exception;
        return undef;
    }
#     print length($aligned_allele), " ", substr($ref_aligned_seq, 0, length($aligned_allele)), "\n";

    ## The get_sub_sorted_alignment does not set the aligned_sequence keys, this is used later on to
    ## mark that these sequences have not been re-aligned yet.
    my $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, $trim_length_in_cols,
            $alignment_length);
#     print_sorted_alignment($sub_sorted_alignment);
    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
    if (@{$sub_sorted_alignment->{sequences}} < 2) {
        ## No sequences left apart from the reference one.
        ## Don't mark it as an excpetion but as missing coverage only
        $EXCEPTION = "REF.SEQ.ONLY";
        return undef;
    }
#     print_sorted_alignment($sub_sorted_alignment);

    return ($sub_sorted_alignment, $alt_allele);
}

my $_last_emf_original_alignment;
sub get_emf_original_alignment {
    my ($this_emf_file, $offset) = @_;
    my $sorted_alignment;

    if ($_last_emf_original_alignment->{$this_emf_file}->{$offset}) {
        return $_last_emf_original_alignment->{$this_emf_file}->{$offset};
    } else {
        $_last_emf_original_alignment = undef;
    }

    $this_emf_file =~ s/\.gz$//;

    my $emf_fh;
    if (-e "$emf/$this_emf_file") {
        open(EMF, "$emf/$this_emf_file") or die "Cannot open file $this_emf_file\n";
        $emf_fh = *EMF;
    } elsif (-e "$emf/$this_emf_file.gz") {
        unless ($allow_gzip) {
            die "Please, use uncomressed EMF files for efficiency\n";
        }
        use IO::Uncompress::Gunzip;
        $emf_fh = new IO::Uncompress::Gunzip "$emf/$this_emf_file.gz" or die "";
    }
    seek($emf_fh, $offset, 0);
    $sorted_alignment = read_sorted_alignment_from_emf_fh($emf_fh);
#     print_sorted_alignment($sorted_alignment);

    $_last_emf_original_alignment->{"$this_emf_file.gz"}->{$offset} = $sorted_alignment;
    close($emf_fh);

    return $sorted_alignment;
}
