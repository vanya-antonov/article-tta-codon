#!/usr/bin/perl -w

use strict;
use warnings;

##################################################################
# Antonov Ivan                  ver.1.00
#

use Getopt::Long;
use Data::Dumper;

##################################################################
# CONSTANTS
my $EXIST   = 0;
my $ID_COLS = undef;

###################################################################
# Parse input data
GetOptions(
	'invert'     => \$EXIST,
	'id_cols=s'  => \$ID_COLS,
) or die usage();

die usage() if @ARGV < 2;

############
my $START_TIME = time;
run(
	all_files  => \@ARGV,
	id_cols    => $ID_COLS ? [ split(/,/, $ID_COLS) ] : [ map { 1 } @ARGV ],
	exist      => $EXIST,
	verbose    => 1,
);
warn "\nElapsed time: ". (time-$START_TIME) ."\n";
############

###################################################################
# SUBROUTINES
sub run
{
	my %opts =@_;
	
	die "The number of files and ID_COLs do not match: ".Dumper(\%opts) if scalar(@{$opts{all_files}}) != scalar(@{$opts{id_cols}});
	
	my $fn_1  = shift @{$opts{all_files}};
	my $col_1 = shift @{$opts{id_cols}};
	
	# Get ids from file2, ...
	my $ids_2 = {};
	while( my $fn = shift @{$opts{all_files}} )
	{
		warn "Reading the file '$fn'...\n" if $opts{verbose};
		my $id_col = shift @{$opts{id_cols}};
		
		open( F , $fn ) or die "Can't open file '$fn': $!";
		while( <F> )
		{
			my @fields = split( /\t/ , $_ );
			die "Wrong column number ($id_col). There are '".scalar(@fields)."' fields only!" if scalar(@fields) < $id_col;
			
			$fields[ $id_col - 1 ] =~ s/^\s*//;
			$fields[ $id_col - 1 ] =~ s/\s*$//;
			
			$ids_2->{ uc($fields[$id_col-1]) } = 1;
		}
		close F;
	}
	
	# Read file_1
	open( F , $fn_1 ) or die "Can't open file '$fn_1': $!";
	while( my $line = <F> )
	{
		next if $line =~ /^\s*$/;
		
		my @fields = split( /\t/ , $line );
		die "Wrong column1 ($col_1). There are '".scalar(@fields)."' fields only!" if scalar(@fields) < $col_1;
		
		my $id_1 = $fields[ $col_1 - 1 ];
		
		$id_1 =~ s/^\s*//;
		$id_1 =~ s/\s*$//;
		
		if( $ids_2->{ uc($id_1) } )
		{
			# $id_1 exists in file2.txt
			print $line if $opts{exist};
		}
		else
		{
			# $id_1 doesn't exist in file2.txt
			print $line unless $opts{exist};
		}
	}
	close F;
	
};

sub usage
{
	my ($script) = ($0 =~ /.*[\/\\](.*)/);
	return "
DESCRIPTION:
     Prints all lines from FILE_1.txt with ids that do NOT appear in
     files FILE_2.txt [FILE_3.txt ... FILE_N.txt]

USAGE:
     $script   [OPTIONS]   FILE_1.txt   FILE_2.txt   [FILE_3.txt ...]   >   OUT.txt

OPTIONS:
     --invert       --  opposite functionality: print all lines from
                        FILE_1.txt with ids that APPEAR in other file(s)
     --id_cols <S>  --  a string with N numbers separated by commas. Default: 1,[...,]1

";
}

