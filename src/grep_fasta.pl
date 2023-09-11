#!/usr/bin/perl -w

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#

$| = 1; # Turn off buffering

use Data::Dumper;
use Getopt::Long;

use MyLib::BaseUtil qw(revcomp);
use MyLib::classes::dataTransform::Fasta_formatManager qw(read_next_seq);

##################################################################
# CONSTANTS
my $SEQ_RE         = '';
my $EXCLUDE_SEQ_RE = '';
my $SEQ_SHORTER    = '';
my $EXACT_MATCH    = 0;
my $INVERT         = 0;
my $PATTERN_STR    = 0;
my $ONE_SEQ_PER_RE = 0;
my $ACGT_FILTER    = 0;

###################################################################
# Parse input data
GetOptions(
	'seq_re=s'            => \$SEQ_RE,
	'exclude_seq_re=s'    => \$EXCLUDE_SEQ_RE,
	'exact_match'         => \$EXACT_MATCH,
	'invert'              => \$INVERT,
	'pattern_str'         => \$PATTERN_STR,
	'seq_shorter=i'       => \$SEQ_SHORTER,
	'one_seq_per_RE'      => \$ONE_SEQ_PER_RE,
	'acgt_filter'         => \$ACGT_FILTER,
) || die usage();

die usage() if @ARGV!=2;

############
my $START_TIME = time;

run(
	patterns       => $ARGV[0],
	fasta_fn       => $ARGV[1],
	seq_re         => $SEQ_RE ? qr/$SEQ_RE/ : '',
	exclude_seq_re => $EXCLUDE_SEQ_RE ? qr/$EXCLUDE_SEQ_RE/ : '',
	exact_match    => $EXACT_MATCH,
	invert         => $INVERT,
	pattern_str    => $PATTERN_STR,
	seq_shorter    => $SEQ_SHORTER,
	one_seq_per_re => $ONE_SEQ_PER_RE,
	ACGT_filter    => $ACGT_FILTER,
);

warn "\nElapsed time: ".(time-$START_TIME)." sec\n";
############

###################################################################
# SUBROUTINES
sub run
{
	my %opts = @_;
	
	my $all_re = {};
	if( $opts{pattern_str} )
	{
		$all_re = { $opts{patterns} => { re => qr/$opts{patterns}/ } };
	}
	else
	{
		$all_re = read_patterns($opts{patterns}, %opts);
	}
	
	# Grep fasta
	open(my $fh, '<', $opts{fasta_fn}) || die "Can't open file $opts{fasta_fn}: $!";
	while(my $seq = read_next_seq($fh))
	{
		print STDERR "\r$seq->{fullname}                 ";
		next if $opts{ACGT_filter}    && $seq->{seq} =~ /[^ACGT]/i;
		next if $opts{exclude_seq_re} && $seq->{seq} =~ $opts{exclude_seq_re};
		next if $opts{seq_re}         && $seq->{seq} !~ $opts{seq_re};
		next if $opts{seq_shorter}    && length($seq->{seq}) >  $opts{seq_shorter};
		last if scalar( keys %$all_re ) == 0;
		
		my $match = 0;
		if( $opts{exact_match} )
		{
			if( $all_re->{ $seq->{fullname} } )
			{
				$match = 1;
				delete $all_re->{ $seq->{fullname} } if $opts{one_seq_per_re};
			}
		}
		else
		{
			foreach my $re_key ( keys %$all_re )
			{
				my $item = $all_re->{$re_key};
				if($seq->{fullname} =~ $item->{re})
				{
					$match = 1;
					delete $all_re->{$re_key} if $opts{one_seq_per_re};
					last;
				}
			}
		}
		
		if( ($match && !$opts{invert}) || (!$match && $opts{invert}) )
		{
			print ">$seq->{fullname}\n$seq->{seq}\n";
		}
	}
	print STDERR "\n";
	close $fh;
}

###
# Returns array of regular expressions
sub read_patterns
{
	my($fn, %opts) = @_;
	
	open(my $fh, $fn eq '-' ? '<&STDIN' : $fn) || die "Can't open file $fn: $!";
    my($pats, $num_cols) = ({}, undef);
	while( my $line = <$fh> )
	{
		$line =~ s/[\n\r]//g;
		next if $line =~ /^\s*$/;
		
		my @vals = split /\t/, $line;
		$num_cols = scalar(@vals) if !defined $num_cols;
		die "Wrong number of columns in the table!" if scalar(@vals) != $num_cols;
		
		if( $opts{exact_match} )
		{
			$pats->{$vals[0]} = 1;
		}
		else
		{
			$vals[0] =~ s/^\s*(.*?)\s*$/$1/;
			$pats->{$vals[0]} = { re => qr/$vals[0]/ };
		}
	}
	close $fh;
	
	return $pats;
}

sub usage
{
	my($script) = $0 =~ /([^\\\/]+)$/;
	return"
DESCRIPTION:
    This is analog of UNIX's grep function that is adapted to files in fasta format.
    Selection of sequences performed by applying patters from <patterns.txt> to each
    info-line in fasta file. If at least one pattern matches the line this sequence
    is selected and will appear in <res.mfa>.
    Patterns are the Perl regular expressions, i.e. each line is used in the following way: /\$LINE/.
    <patterns.txt> is a single-column table without header.

USAGE:
    $script   <patterns.txt>   <seqs.mfa>   >   <res.mfa>

OPTIONS:
    --pattern_str               --  <patterns.txt> is a single pattern string rather than file with patterns
    --seq_re <reg_exp>          --  regular expression to apply to sequence (not sequence name). The sequence
                                    will be printed out only if sequence name matches <patterns.txt>
                                    and the sequence itsefl matches <seq_re>.
    --exclude_seq_re <reg_exp>  --  regular expression to apply to sequence (not sequence name).
                                    If this RE matches a sequence, the sequence will not be
                                    printed (even if there were matched pattern from <patterns.txt>).
    --exact_match               --  exact match is required (use simple string matching intead of RE)
    --invert                    --  Invert the sense of matching, to select non-matching lines.
                                    So we have opposite functionality -- print all sequences from <seqs.mfa>
                                    EXCEPT those that match at least one RE from <patterns.txt>
    --seq_shorter <NUM>         --  output sequences shorter than <NUM> only
    --one_seq_per_RE            --  optimization option. The reg-expr is removed after the first match.
    --ACGT_filter               --  do not output sequences with non-standard letter even if they match a pattern

";
}

