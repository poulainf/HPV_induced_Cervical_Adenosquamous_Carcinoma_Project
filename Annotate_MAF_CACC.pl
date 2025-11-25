#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
#use Smart::Comments;
use warnings FATAL => 'uninitialized';

# Récupération des arguments
my $infile     = shift;
my $list_exon  = shift;
my $suffix = shift // '';  # suffixe optionnel

if ($suffix eq "INDELs"){
# Lecture des fichiers d’info (référence à un hash)
my $exon_ref  = read_info2($list_exon);

# Définition du fichier de sortie
(my $outname = $infile) =~ s/\.txt$//;   
my $outfile_tot = $outname . "_cadd_ANNOTED" . $suffix . ".txt";

open my $in,  '<', $infile;
open my $out, '>', $outfile_tot;

LINE:
while (my $line = <$in>) {
    chomp $line;
    my @infos = split /\t/, $line;

    if ($line =~ m/Type_score/xms){
        say {$out} join "\t",$line,"CADD_PHRED_1.7";
        next LINE;
    }

    my $chr      = $infos[16];
    my $position = $infos[187];
    my $Ref      = $infos[189];
    my $Alt      = $infos[190];

    $chr =~ s/^chr//;  
    

	if ( $Ref eq '-' ) {
		($Ref) = $Alt =~ /^(.).*/;
	}
	
	my $len_Alt = length $Alt;

    #my $seg_ID = join "_", $chr, $position, $Ref, $len_Alt;
	my $seg_ID = join "_", $chr, $position, $Ref, $Alt;


    if (defined (my $Phred = $exon_ref->{$seg_ID})) {
        ## $Phred
        
            say {$out} join "\t", $line, $Phred;		  
        
    }
}
close $in;
close $out;
}else{
# Lecture des fichiers d’info (référence à un hash)
my $exon_ref  = read_info2($list_exon);

# Définition du fichier de sortie
(my $outname = $infile) =~ s/\.txt$//;   
my $outfile_tot = $outname . "_cadd_ANNOTED" . $suffix . ".txt";

open my $in,  '<', $infile;
open my $out, '>', $outfile_tot;

LINE:
while (my $line = <$in>) {
    chomp $line;
    my @infos = split /\t/, $line;

    if ($line =~ m/Type_score/xms){
        say {$out} join "\t",$line,"CADD_PHRED_1.7";
        next LINE;
    }

    my $chr      = $infos[16];
    my $position = $infos[17];
    my $Ref      = $infos[21];
    my $Alt      = $infos[22];

    $chr =~ s/^chr//;  
    

	if ( $Ref eq '-' ) {
		($Ref) = $Alt =~ /^(.).*/;
	}
	if ( $Alt eq '-' ) {
		($Alt) = $Ref =~ /^(.).*/;
	}

    my $seg_ID = join "_", $chr, $position, $Ref, $Alt;

    if (defined (my $Phred = $exon_ref->{$seg_ID})) {
        ## $Phred
        
            say {$out} join "\t", $line, $Phred;		  
        
    }
}
close $in;
close $out;
}

###############################  

sub read_info2 {
    my $infile3 = shift;
    my %Read_info2;

    open my $fh, '<', $infile3;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ m/^#/;  # ignorer lignes de commentaire

        my @infos3 = split /\t/, $line;
        my ($chr, $position, $Ref, $Alt, $RawScore, $PHRED) = @infos3[0..5];

        my $id = join "_", $chr, $position, $Ref, $Alt;
        ## $id
        $Read_info2{$id} = $PHRED;
    }
    return \%Read_info2;   # retour = référence
}


sub read_info3 {
    my $infile4 = shift;
    my %Read_info3;

    open my $fh, '<', $infile4;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ m/^#/;  # ignorer lignes de commentaire

        my @infos3 = split /\t/, $line;
        my ($chr, $position, $Ref, $Alt, $RawScore, $PHRED) = @infos3[0..5];
		my $len_Alt = length $Alt;
        my $id = join "_", $chr, $position, $Ref, $len_Alt;
        ## $id
        $Read_info3{$id} = $PHRED;
    }
    return \%Read_info3;   # retour = référence
}
