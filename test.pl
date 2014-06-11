#! /usr/bin/env perl
use strict;
use warnings;

use Test::More;

use BaseAncestralAlleles;
use RunAncestralAllelesOnEMF;

my $normal_flank = "ACTGACTGAC";
my $all_N_flank = "NNNNNNNNNN";
my $low_complexity_flank = "ACACACACAC";

subtest "get_reference_sequence_exception" => sub {
    
    ## get_reference_sequence_excpetion($flank5, $flank3, $alignment_length, $max_alignment_length);
    is(get_reference_sequence_exception($normal_flank, $normal_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($normal_flank, $normal_flank, 101, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($all_N_flank, $normal_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($normal_flank, $all_N_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($all_N_flank, $all_N_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($low_complexity_flank, $normal_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($normal_flank, $low_complexity_flank, 10, 100, "/tmp/"), undef);
    isnt(get_reference_sequence_exception($low_complexity_flank, $low_complexity_flank, 10, 100, "/tmp/"), undef);
};


my $seq1 = {
    aligned_sequence => "ACTGACTGACTG",
    species => "s1",
    chr => "c1",
    start => 121,
    end => 132,
    strand => 1,
};
my $seq2 = {
    aligned_sequence => "ACTGACTGACTG",
    species => "s2",
    chr => "c2",
    start => 221,
    end => 232,
    strand => 1,
};
my $seq3 = {
    aligned_sequence => "ACTGACTGACTG",
    species => "s3",
    chr => "c3",
    start => 321,
    end => 332,
    strand => 1,
};
my $seq4 = {
    aligned_sequence => "ACTGACTGACTG",
    species => "s4",
    chr => "c4",
    start => 421,
    end => 432,
    strand => 1,
};

my $sorted_alignment = {
    tree => '(((1:0.1,0:0.2):0.3,2:0.4):0.5,3:0.6):0.7;',
    positions => [1, 0, 2, 3],
    sequences => [$seq1, $seq2, $seq3, $seq4],
};

my $sub_sorted_alignment;

## Easy sub-alignment

subtest "get_sub_sorted_alignment (easy)" => sub {
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 0, 10);
    
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "ACTGACTGAC");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 121);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 130);
};

## Reverse strand

subtest "get_sub_sorted_alignment (reverse strand)" => sub {
    
    $sorted_alignment->{sequences}->[0]->{strand} = -1;
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 0, 10);
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "ACTGACTGAC");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 123);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 132);
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 1, 10);
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "CTGACTGACT");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 122);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 131);
};

## Gapped alignment

subtest "get_sub_sorted_alignment (gapped alignment)" => sub {
    
    $sorted_alignment->{sequences}->[0]->{strand} = 1;
    $sorted_alignment->{sequences}->[0]->{aligned_sequence} = "A----CTGACTG";
    $sorted_alignment->{sequences}->[0]->{end} = 128;
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 0, 10);
    
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "ACTGAC");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 121);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 126);
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 1, 10);
    
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "CTGACT");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 122);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 127);
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "CTGACTG");
    is($sub_sorted_alignment->{sequences}->[0]->{start}, 122);
    is($sub_sorted_alignment->{sequences}->[0]->{end}, 128);
};


#subtest "Splicing first sequence" => sub {
#    
#    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
#    $sub_sorted_alignment->{sequences}->[0]->{original_sequence} = "";
#    
#    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
#    
#    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s2");
#    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s3");
#    is($sub_sorted_alignment->{sequences}->[2]->{species}, "s4");
#    is($sub_sorted_alignment->{tree}, "((0:0.5,2:0.4):0.5,3:0.6):0.7;");
#    is(join("-", @{$sub_sorted_alignment->{positions}}), "0-1-2");
#};

subtest "Splicing second sequence" => sub {
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    $sub_sorted_alignment->{sequences}->[1]->{original_sequence} = "";
    
    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
    
    is(@{$sub_sorted_alignment->{sequences}}, 3);
    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s1");
    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s3");
    is($sub_sorted_alignment->{sequences}->[2]->{species}, "s4");
    is($sub_sorted_alignment->{tree}, "((0:0.5,1:0.4):0.5,2:0.6):0.7;");
    is(join("-", @{$sub_sorted_alignment->{positions}}), "0-1-2");
};

subtest "Splicing third sequence" => sub {
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    $sub_sorted_alignment->{sequences}->[2]->{original_sequence} = "";
    
    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
    
    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s1");
    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s2");
    is($sub_sorted_alignment->{sequences}->[2]->{species}, "s4");
    is($sub_sorted_alignment->{tree}, "((1:0.1,0:0.2):0.8,2:0.6):0.7;");
    is(join("-", @{$sub_sorted_alignment->{positions}}), "1-0-2");
};

subtest "Splicing fourth sequence" => sub {
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    $sub_sorted_alignment->{sequences}->[3]->{original_sequence} = "";
    
    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
    
    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s1");
    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s2");
    is($sub_sorted_alignment->{sequences}->[2]->{species}, "s3");
    is($sub_sorted_alignment->{tree}, "((1:0.1,0:0.2):0.3,2:0.4):1.2;");
    is(join("-", @{$sub_sorted_alignment->{positions}}), "1-0-2");
};

#subtest "Splicing first and second sequences" => sub {
#    
#    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
#    $sub_sorted_alignment->{sequences}->[0]->{original_sequence} = "";
#    $sub_sorted_alignment->{sequences}->[1]->{original_sequence} = "";
#    
#    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
#    
#    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s3");
#    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s4");
#    is($sub_sorted_alignment->{tree}, "(2:0.9,3:0.6):0.7;");
#    is(join("-", @{$sub_sorted_alignment->{positions}}), "0-1");
#};

subtest "Splicing third and fourth sequences" => sub {

    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    $sub_sorted_alignment->{sequences}->[2]->{original_sequence} = "";
    $sub_sorted_alignment->{sequences}->[3]->{original_sequence} = "";
    
    splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
    
    is($sub_sorted_alignment->{sequences}->[0]->{species}, "s1");
    is($sub_sorted_alignment->{sequences}->[1]->{species}, "s2");
    is($sub_sorted_alignment->{tree}, "(1:0.1,0:0.2):1.5;");
    is(join("-", @{$sub_sorted_alignment->{positions}}), "1-0");
};


subtest 'reverse complement sorted_alignment' => sub {

    $sorted_alignment->{sequences}->[0]->{strand} = -1;
    
    straighten_sorted_alignment($sorted_alignment);
    
    is($sorted_alignment->{sequences}->[0]->{strand}, 1);
    is($sorted_alignment->{sequences}->[1]->{strand}, -1);
    is($sorted_alignment->{sequences}->[2]->{strand}, -1);
    is($sorted_alignment->{sequences}->[3]->{strand}, -1);
    is($sorted_alignment->{sequences}->[0]->{aligned_sequence}, "CAGTCAG----T");
    is($sorted_alignment->{sequences}->[1]->{aligned_sequence}, "CAGTCAGTCAGT");
    is($sorted_alignment->{sequences}->[2]->{aligned_sequence}, "CAGTCAGTCAGT");
    is($sorted_alignment->{sequences}->[3]->{aligned_sequence}, "CAGTCAGTCAGT");
    
    $sorted_alignment->{sequences}->[0]->{strand} = -1;
    straighten_sorted_alignment($sorted_alignment); # undo
};

subtest 'reverse complement sub_sorted_alignment' => sub {
    
    $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, 2, 10);
    $sub_sorted_alignment->{sequences}->[0]->{strand} = -1;
    
    straighten_sorted_alignment($sub_sorted_alignment);
    
    is($sub_sorted_alignment->{sequences}->[0]->{strand}, 1);
    is($sub_sorted_alignment->{sequences}->[1]->{strand}, -1);
    is($sub_sorted_alignment->{sequences}->[2]->{strand}, -1);
    is($sub_sorted_alignment->{sequences}->[3]->{strand}, -1);
    is($sub_sorted_alignment->{sequences}->[0]->{original_sequence}, "CAGTCAG");
    is($sub_sorted_alignment->{sequences}->[1]->{original_sequence}, "CAGTCAGTCA");
    is($sub_sorted_alignment->{sequences}->[2]->{original_sequence}, "CAGTCAGTCA");
    is($sub_sorted_alignment->{sequences}->[3]->{original_sequence}, "CAGTCAGTCA");
};

my $ortheus_alignment;
my ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3);

subtest 'Alleles from easy deletion (ref seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_deletion($ortheus_alignment, 10, "A");
    
    is($ref_allele, "A");
    is($sis_allele, "A");
    is($anc_allele, "A");
    is($old_allele, "A");
    is($length_flank5, 10);
    is($length_allele, 1);
    is($length_flank3, 9);
};

subtest 'Alleles from easy deletion (alt seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGAC-CTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_deletion($ortheus_alignment, 10, "A");
    
    is($ref_allele, "");
    is($sis_allele, "A");
    is($anc_allele, "A");
    is($old_allele, "A");
    is($length_flank5, 10);
    is($length_allele, 1);
    is($length_flank3, 9);
};

subtest 'Alleles from deletion on 3bp homopolymer (ref seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_deletion($ortheus_alignment, 10, "A");
    
    is($ref_allele, "AAA");
    is($sis_allele, "AAA");
    is($anc_allele, "AAA");
    is($old_allele, "AAA");
    is($length_flank5, 8);
    is($length_allele, 3);
    is($length_flank3, 9);
};

subtest 'Alleles from deletion on 2bp homopolymer (alt seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGAA-CTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGAAACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_deletion($ortheus_alignment, 10, "A");
    
    is($ref_allele, "AA");
    is($sis_allele, "AAA");
    is($anc_allele, "AAA");
    is($old_allele, "AAA");
    is($length_flank5, 8);
    is($length_allele, 3);
    is($length_flank3, 9);
};

subtest 'Alleles from easy insertion (ref seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_insertion($ortheus_alignment, 10, "G");
    
    is($ref_allele, "");
    is($sis_allele, "");
    is($anc_allele, "");
    is($old_allele, "");
    is($length_flank5, 10);
    is($length_allele, 0);
    is($length_flank3, 10);
};


subtest 'Alleles from easy insertion (alt seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGACGACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_insertion($ortheus_alignment, 10, "G");
    
    is($ref_allele, "G");
    is($sis_allele, "");
    is($anc_allele, "");
    is($old_allele, "");
    is($length_flank5, 10);
    is($length_allele, 1);
    is($length_flank3, 10);
};

subtest 'Alleles from insertion next to same nucleotide (ref seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGACACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_insertion($ortheus_alignment, 10, "C");
    
    is($ref_allele, "C");
    is($sis_allele, "C");
    is($anc_allele, "C");
    is($old_allele, "C");
    is($length_flank5, 9);
    is($length_allele, 1);
    is($length_flank3, 10);
};

subtest 'Alleles from insertion next to same nucleotide (alt seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "ACTGACTGACCACTGACTGAC" },
        'sis' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
        'anc' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
        'old' => { 'aligned_sequence' => "ACTGACTGAC-ACTGACTGAC" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_insertion($ortheus_alignment, 10, "C");
    
    is($ref_allele, "CC");
    is($sis_allele, "C");
    is($anc_allele, "C");
    is($old_allele, "C");
    is($length_flank5, 9);
    is($length_allele, 2);
    is($length_flank3, 10);
};


subtest 'Alleles from insertion next to same nucleotide (alt seq)' => sub {
    
    $ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "--------------TGGATAAGTGTAGGCAATTT" },
        'sis' => { 'aligned_sequence' => "--------------TGGATAAGTGTAGGCAATTT" },
        'anc' => { 'aligned_sequence' => "--------------TGGATAAGTGTAGGCAATTT" },
        'old' => { 'aligned_sequence' => "--------------TGGATAAGTGTAGGCAATTT" },
    };
    
    ($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3) = BaseAncestralAlleles::get_alleles_for_insertion($ortheus_alignment, 10, "C");
    
    is($ref_allele, "");
    is($sis_allele, "");
    is($anc_allele, "");
    is($old_allele, "");
    is($length_flank5, 24);
    is($length_allele, 0);
    is($length_flank3, 10);
};


subtest 'call_ancestral_allele_for_deletion -> deletion' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
    };
    
    my $deletion_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGG-GCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
    };
    
    my ($ref_reference_allele, $ref_deletion_allele, $ancestral_allele_call, $indel_call) = call_ancestral_allele_for_deletion(undef, $reference_ortheus_alignment, $deletion_ortheus_alignment, "G", 10, 0);
    
    is($ref_reference_allele, "GGGGG");
    is($ref_deletion_allele, "GGGG");
    is($ancestral_allele_call, "GGGgG");
    is($indel_call, "deletion");
};

subtest 'call_ancestral_allele_for_deletion -> insertion' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGG-TGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGG-GGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGG-GGCAGTG" },
    };
    
    my $deletion_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGGCAGTG" },
    };
    
    my ($ref_reference_allele, $ref_deletion_allele, $ancestral_allele_call, $indel_call) = call_ancestral_allele_for_deletion(undef, $reference_ortheus_alignment, $deletion_ortheus_alignment, "G", 10, 0);
    
    is($ref_reference_allele, "GGGGG");
    is($ref_deletion_allele, "GGGG");
    is($ancestral_allele_call, "GGgG");
    is($indel_call, "insertion");
};


subtest 'call_ancestral_allele_for_deletion -> complex_deletion' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
    };
    
    my $deletion_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGG-GCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTG" },
    };
    
    my ($ref_reference_allele, $ref_deletion_allele, $ancestral_allele_call, $indel_call) = call_ancestral_allele_for_deletion(undef, $reference_ortheus_alignment, $deletion_ortheus_alignment, "G", 10, 0);
    
    is($ref_reference_allele, "GGGGG");
    is($ref_deletion_allele, "GGGG");
    is($ancestral_allele_call, "GGGTG");
    is($indel_call, "complex_deletion");
};

subtest 'call_ancestral_allele_for_deletion -> complex_insertion' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGG-TGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGG-TGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGG-TGCAGTG" },
    };
    
    my $deletion_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGCAGTG" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGTGCAGTG" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGTGCAGTG" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGTGCAGTG" },
    };
    
    my ($ref_reference_allele, $ref_deletion_allele, $ancestral_allele_call, $indel_call) = call_ancestral_allele_for_deletion(undef, $reference_ortheus_alignment, $deletion_ortheus_alignment, "G", 10, 0);
    
    is($ref_reference_allele, "GGGGG");
    is($ref_deletion_allele, "GGGG");
    is($ancestral_allele_call, "GGTG");
    is($indel_call, "complex_insertion");
};

subtest 'call_ancestral_allele_for_insertion' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGG-GCAGTGA" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTGA" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTGA" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTGA" },
    };
    
    my $insertion_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTGA" },
        'sis' => { 'aligned_sequence' => "CATTAAGGCTGGGTGCAGTGA" },
        'anc' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTGA" },
        'old' => { 'aligned_sequence' => "CATTAAGGCTGGGGGCAGTGA" },
    };
    
    my ($ref_reference_allele, $ref_insertion_allele, $ancestral_allele_call, $indel_call) = call_ancestral_allele_for_insertion(undef, $reference_ortheus_alignment, $insertion_ortheus_alignment, "G", 10, 0);
    
    is($ref_reference_allele, "GGGG");
    is($ref_insertion_allele, "GGGGG");
    is($ancestral_allele_call, "GGGgG");
    is($indel_call, "deletion");
};

subtest 'complex splice tree' => sub {
    my $sorted_alignment = {
        'tree' => '((((((4:0.0066,5:0.0084):0,3:0.0098):0.0010,6:0.0225):0.0010,7:0.0429):0.0010,((0:0.0049,1:0.0060):0.0076,2:0.0123):0.0309):0.0448,8:0.0438):0.0010;',
        'positions' => [
        '4',
        '5',
        '3',
        '6',
        '7',
        '0',
        '1',
        '2',
        '8'
        ],
        'sequences' => [
        {
            'original_sequence' => 'TCGTCTCTACTAAAAAAACA',
            'chr' => '22',
            'name' => 'Hsap_22_25795554_25795573[+]',
            'strand' => 1,
            'species' => 'homo_sapiens',
            'end' => 25795573,
            'start' => 25795554
        },
        {
            'original_sequence' => 'TCGTCTCTACTAAAAAAACA',
            'chr' => '22',
            'name' => 'Ptro_22_23945479_23945498[+]',
            'strand' => 1,
            'species' => 'pan_troglodytes',
            'end' => 23945498,
            'start' => 23945479
        },
        {
            'original_sequence' => 'TCGTCTCTACTAAAAAAACA',
            'chr' => '22',
            'name' => 'Ggor_22_9399369_9399388[+]',
            'strand' => 1,
            'species' => 'gorilla_gorilla',
            'end' => 9399388,
            'start' => 9399369
        },
        {
            'original_sequence' => '',
            'chr' => '11',
            'name' => 'Ggor_11_65264918_65264917[-]',
            'strand' => -1,
            'species' => 'gorilla_gorilla',
            'start' => 65264918,
            'end' => 65264917
        },
        {
            'original_sequence' => '',
            'chr' => '11',
            'name' => 'Hsap_11_68085445_68085444[-]',
            'strand' => -1,
            'species' => 'homo_sapiens',
            'start' => 68085445,
            'end' => 68085444
        },
        {
            'original_sequence' => '',
            'chr' => '11',
            'name' => 'Ptro_11_66013089_66013088[-]',
            'strand' => -1,
            'species' => 'pan_troglodytes',
            'start' => 66013089,
            'end' => 66013088
        },
        {
            'original_sequence' => '',
            'chr' => '11',
            'name' => 'Pabe_11_7800984_7800983[+]',
            'strand' => 1,
            'species' => 'pongo_abelii',
            'end' => 7800983,
            'start' => 7800984
        },
        {
            'original_sequence' => '',
            'chr' => '14',
            'name' => 'Mmul_14_6398672_6398671[+]',
            'strand' => 1,
            'species' => 'macaca_mulatta',
            'end' => 6398671,
            'start' => 6398672
        },
        {
            'original_sequence' => '',
            'chr' => '11',
            'name' => 'Cjac_11_123845304_123845303[-]',
            'strand' => -1,
            'species' => 'callithrix_jacchus',
            'start' => 123845304,
            'end' => 123845303
        }
        ]
    };
    
    splice_uninformative_sequences_from_alignment($sorted_alignment);

    is($sorted_alignment->{tree}, '((0:0.0049,1:0.0060):0.0076,2:0.0123):0.0767;');
};

subtest 'complex splice tree' => sub {
    my $sorted_alignment = {
        'tree' => '((4:0,5:0.0001):0.0268,(((2:0.0111,1:0.0114):0,0:0.0114):0.0010,3:0.0263):0.0278):0.0010;',
        'positions' => [
        '4',
        '5',
        '2',
        '1',
        '0',
        '3'
        ],
        'sequences' => [
        {
            'original_sequence' => 'CAACCCCTCATGGACCAGAC',
            'chr' => '22',
            'name' => 'Hsap_22_48652778_48652797[+]',
            'strand' => 1,
            'species' => 'homo_sapiens',
            'end' => 48652797,
            'start' => 48652778
        },
        {
            'original_sequence' => '',
            'chr' => '22',
            'name' => 'Ggor_22_32825866_32825865[+]',
            'strand' => 1,
            'species' => 'gorilla_gorilla',
            'end' => 32825865,
            'start' => 32825866
        },
        {
            'original_sequence' => 'CAACCCCTCATGGACCAGAC',
            'chr' => '22',
            'name' => 'Ptro_22_47085124_47085143[+]',
            'strand' => 1,
            'species' => 'pan_troglodytes',
            'end' => 47085143,
            'start' => 47085124
        },
        {
            'original_sequence' => 'CAACCCCTCATGGACCAG',
            'chr' => '22',
            'name' => 'Pabe_22_43745290_43745307[+]',
            'strand' => 1,
            'species' => 'pongo_abelii',
            'end' => 43745307,
            'start' => 43745290
        },
        {
            'original_sequence' => '',
            'chr' => '10',
            'name' => 'Mmul_10_92272637_92272636[+]',
            'strand' => 1,
            'species' => 'macaca_mulatta',
            'end' => 92272636,
            'start' => 92272637
        },
        {
            'original_sequence' => '',
            'chr' => '1099214733268',
            'name' => 'Mmul_1099214733268_6628_6627[+]',
            'strand' => 1,
            'species' => 'macaca_mulatta',
            'end' => 6627,
            'start' => 6628
        }
        ]
    };

    splice_uninformative_sequences_from_alignment($sorted_alignment);
    
    is($sorted_alignment->{tree}, '((1:0.0111,0:0.0114):0.0010,2:0.0263):0.0288;');
};


subtest 'Bug in extracting the flanks (1)' => sub {
    
    my $reference_ortheus_alignment = {
        'ref' => { 'aligned_sequence' => "TGAAGTGCCT-GGTCAGCTTA" },
        'sis' => { 'aligned_sequence' => "TGAAGTGCCT-GGTCAGCTTA" },
        'anc' => { 'aligned_sequence' => "TGAAGTGCCT-GGTCAGCTTA" },
        'old' => { 'aligned_sequence' => "TGAAGTGCCT-GGTCAGCTTA" },
    };
    BaseAncestralAlleles::get_alleles_for_insertion($reference_ortheus_alignment, 10, "A");
    BaseAncestralAlleles::get_alleles_for_insertion($reference_ortheus_alignment, 10, "C");
    BaseAncestralAlleles::get_alleles_for_insertion($reference_ortheus_alignment, 10, "T");
    BaseAncestralAlleles::get_alleles_for_insertion($reference_ortheus_alignment, 10, "G");
    pass();
};


#GCA---CCCACAA-CATGCCCAG

done_testing();
