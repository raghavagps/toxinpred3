#!/usr/bin/perl -w

###############################################################################################
# This is the MERCI (Motif - EmeRging and with Classes - Identification) motif locator program.
# Copyright (C) 2010  Celine Vens
# 
# Usage will be shown when script is run without any arguments: "perl MERCI_motif_locator.pl"
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#     
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#     
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################################    


use strict;
use Time::HiRes;


# GLOBAL VARIABLES
##################

# input arguments
my @arguments = @ARGV;
my $posfile = -1; # fasta file with positive sequences (REQUIRED)
my $negfile = -1; # fasta file with negative sequences
my $motiffile = "motifs"; # output file to store motifs
my $hierarchy = "NONE"; # class hierarchy to use (NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL/hierarchy_file)
my $stringlength = 10000;
my $outputfile = $motiffile . ".occurrences";
my $maxgaplength = 1;

my %string_representation; # the representation of the concepts in the data file
my %coverage_hash = ();


# PROCESSING ARGUMENTS
######################
#print "\nThis is the MERCI (Motif - EmeRging and with Classes - Identification) motif locator program.\n\n";

if(@arguments <= 0 ){	
	print "To run the program, please define at least the required options:\n\n";
	print "  -p posfile          # the file to search for occurrences of the motifs (REQUIRED)\n";
	print "  -n negfile          # optional second file to search for occurrences of the motifs (e.g. the negative sequences)\n";
	print "  -i inputfile        # the file with motifs (default: motifs)\n";
	print "  -o outputfile       # the file where motif occurrences will be written (default: inputfile.occurrences)\n";
	print "  -c classification   # the classification hierarchy to use (values: NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL/my_classification_file, default: NONE)\n";
	print "  -gl maxgaplength    # maximal gap length (default: 1)\n";
	print "  -s sequencelength   # only the first sequencelength positions of the sequences will be considered (default: 10000)\n\n";
	print "See manual for more information on the options.\n\n";
	exit;
}
else
{
	my $option;
	my $value;
	while ($option = shift(@arguments))
	{
		$value = shift(@arguments);
		&process_option($option,$value);
	}
}
&check_options();

sub process_option()
{
	my $option = $_[0];
	my $value = $_[1];
	
	if ($option eq "-p")
	{
		if (not defined $value)
		{
			die "File $value (option -p) not found.\n";
		}
		elsif (not -e $value)
		{
			die "File $value (option -p) not found.\n";
		}
		$posfile = $value;
	}
	elsif ($option eq "-n")
	{
		if (not defined $value)
		{
			die "File $value (option -n) not found.\n";
		}
		elsif (not -e $value)
		{
			die "File $value (option -n) not found.\n";
		}
		$negfile = $value;
	}
	elsif ($option eq "-i")
	{
		if (not defined $value)
		{
			die "File $value (option -i) not found.\n";
		}
		elsif (not -e $value)
		{
			die "File $value (option -i) not found.\n";
		}
		$motiffile = $value;
		$outputfile = $motiffile . ".occurrences";
	}
	elsif ($option eq "-o")
	{
		$outputfile = $value;
	}
	elsif ($option eq "-c")
	{
		if (((((not $value eq "NONE") and (not $value eq "KOOLMAN-ROHM")) and (not $value eq "BETTS-RUSSELL")) and (not $value eq "RASMOL")))
		{
			if (not -e $value)
			{
				die "Classification file $value (option -c) not found. Use an existing file name, or one of NONE/KOOLMAN-ROHM/BETTS-RUSSELL/RASMOL.\n";
			}
		}
		$hierarchy = $value;
	}
	elsif ($option eq "-s")
	{
		$stringlength = $value;
	}
	elsif ($option eq "-gl")
	{
		$maxgaplength = $value;
	}
	else
	{
		die "Unknown option $option. Run the program without any options to see its usage.\n";
	}
}

sub check_options()
{
	if ($posfile eq "-1")
	{
		die "No fasta file with (positive) sequences given. Please set option -p.\n";
	}
}


# HIERARCHY
###########
#print "   * Loading hierarchy $hierarchy\n";

$string_representation{"root"} = 'root';
$string_representation{"gap"} = '.{0,' . $maxgaplength . '}';

$string_representation{"A"} = 'A';
$string_representation{"C"} = 'C';
$string_representation{"D"} = 'D';
$string_representation{"E"} = 'E';
$string_representation{"F"} = 'F';
$string_representation{"G"} = 'G';
$string_representation{"H"} = 'H';
$string_representation{"I"} = 'I';
$string_representation{"K"} = 'K';
$string_representation{"L"} = 'L';
$string_representation{"M"} = 'M';
$string_representation{"N"} = 'N';
$string_representation{"P"} = 'P';
$string_representation{"Q"} = 'Q';
$string_representation{"R"} = 'R';
$string_representation{"S"} = 'S';
$string_representation{"T"} = 'T';
$string_representation{"U"} = 'U';
$string_representation{"V"} = 'V';
$string_representation{"W"} = 'W';
$string_representation{"Y"} = 'Y';

&load_hierarchy($hierarchy);

sub load_hierarchy()
{
	my $hierarchy = $_[0];
	
	if ($hierarchy eq "NONE") # hierarchy without concepts
	{
		1;
	}
	
	elsif ($hierarchy eq "KOOLMAN-ROHM")
	{
		$string_representation{"aliphatic"} = '[AGILV]';
		$string_representation{"sulfur"} = '[CM]';
		$string_representation{"aromatic"} = '[FYW]';
		$string_representation{"neutral"} = '[STNQ]';
		$string_representation{"acidic"} = '[DE]';
		$string_representation{"basic"} = '[RHK]';
	}
	
	elsif ($hierarchy eq "BETTS-RUSSELL")
	{
		$string_representation{"polar"} = '[HKRDEYWTCSNQ]';
		$string_representation{"charged"} = '[DERHK]';
		$string_representation{"negative"} = '[DE]';
		$string_representation{"positive"} = '[RHK]';
		$string_representation{"small"} = '[AGCSPNDTV]';
		$string_representation{"tiny"} = '[AGCS]';
		$string_representation{"hydrophobic"} = '[HFWYILVMKTAGC]';
		$string_representation{"aromatic"} = '[HFWY]';
		$string_representation{"aliphatic"} = '[ILV]';
	}
	
	elsif ($hierarchy eq "RASMOL")
	{
		$string_representation{"charged"} = '[DERHK]';
		$string_representation{"acidic"} = '[DE]';
		$string_representation{"basic"} = '[RHK]';
		$string_representation{"neutral"} = '[ANCQGILMFPSTWYV]';
		$string_representation{"cyclic"} = '[HFWYP]';
		$string_representation{"acyclic"} = '[ARNDCEQGILKMSTV]';
		$string_representation{"aromatic"} = '[HFWY]';
		$string_representation{"aliphatic"} = '[AGILV]';
		$string_representation{"surface"} = '[RNDEQGHKPSTY]';
		$string_representation{"buried"} = '[ACILMFWV]';
		$string_representation{"hydrophobic"} = '[AGILMFPWYV]';
		$string_representation{"polar"} = '[RNDCEQHKST]';
		$string_representation{"small"} = '[AGS]';
		$string_representation{"medium"} = '[NDCPTV]';
		$string_representation{"large"} = '[REQHILKMFWY]';
	}
	
	else
	{
		&load_user_hierarchy($hierarchy);
	}
}


# main program
##############
# read sequences and store them as an array of strings and headers
my $start = [ Time::HiRes::gettimeofday( ) ];
#print "   * Reading data\n";
my ($posset,$posheaders) = &extract_sequences($posfile); # (global var)
my ($negset,$negheaders);
my $negpresent = 0;
if (not $negfile eq "-1")
{
	$negpresent = 1;
	($negset,$negheaders) = &extract_sequences($negfile); # (global var)
}


open(IN,"$motiffile") or die "Can not open motif input file: $motiffile\n";
open(OUT,">$outputfile");

my $line=<IN>;
while ($line !~ /Motifs:/)
{
	$line=<IN>;
}

while ($line=<IN>)
{
	chomp($line);
	my $pattern = &pattern_to_string($line);
	&check_locations($line,$pattern);
}
close(IN);
close(OUT);

&check_coverage();

#print "   * Motif locations written to $outputfile\n";


# read the input file and return the sequences as an array of strings
sub extract_sequences()
{
	my $file = $_[0];
	open(IN,"$file") or die "file $file not found\n";
	my @lines = <IN>;
	close(IN);
	my @sequenceset;
	my @sequenceheaders;
	my $sequence = "";
	foreach my $l (@lines)
	{
		if ($l =~ /^>/) # header line
		{
			if (not $sequence eq "") # if not first line
			{
				my $substring = substr($sequence,0,$stringlength);
				push (@sequenceset, $substring);
			}
			$sequence = "";
			push (@sequenceheaders,$l);
		}
		else # sequence line, no header line
		{
			chomp($l);
			$l =~ s/^\s+//; # remove whitespace in front and after sequence
			$l =~ s/\s+$//;
			$sequence = $sequence . $l;
		}
	}
	my $substring = substr($sequence,0,$stringlength);
	push (@sequenceset, $substring); # also add last line
	#print "       input file $file read -> " . @sequenceset . " sequences found\n";
	return (\@sequenceset,\@sequenceheaders);
}

# changes a pattern into a string, with the hierarchy concepts replaced by a regular expression
sub pattern_to_string()
{
	my $pattern = $_[0]; # string
	my $patternstring = "";
	my @pattels = split(" ",$pattern);
	foreach my $el (@pattels) # for each ref to string in the array
	{
		if (($el eq "gap") and ($maxgaplength == 0))
		{
			die "Motif found with gap, but maximal gap length is set to zero\n";
		}
		$patternstring = $patternstring . $string_representation{$el};
	}
	return $patternstring;
}

# check at which locations in the sequences the patterns occur
sub check_locations()
{
	my $patternstring = $_[0];
	my $pattern = $_[1];

	my @elements = split(" ", $patternstring);
	my $patternlength = @elements;
	
	print OUT "MOTIF: $patternstring\n";
	print OUT "******\n";
	
	my $count = 0;
	my $containsgap = &contains_gap($patternstring);
	for (my $i=0; $i<@$posset; $i++)
	{
		my $target = ${$posset}[$i];
		my @positions = ();
		my @strings = ();
		if ($containsgap < 1)
		{
			while ($target =~ /$pattern/g)
			{
				my $nextpos = pos $target; # pos points to the next character
				$nextpos = $nextpos-$patternlength+1; # +1 to start counting at 1 instead of 0
				push(@positions,$nextpos);
				my $motif = substr($target,$nextpos,$patternlength);
				push(@strings,$motif);
				pos $target = $nextpos + 1;
			}
			if (@positions > 0)
			{
				my $seqid = $i+1; # to start counting at 1 instead of 0
				print OUT "${$posheaders}[$i] (sequence $seqid)\n";
				print OUT "    position(s): @positions\n";
				print OUT "    motif(s): @strings\n";
				$count++;
			}
		}
		else
		{
			my $found = 0;
			if ($target =~ /$pattern/)
			{
				$found=1;
			}
			if ($found > 0)
			{
				my $seqid = $i+1; # to start counting at 1 instead of 0
				print OUT "${$posheaders}[$i] (sequence $seqid)\n";
				$count++;
			}
		}
		
	}
	print OUT "$count positive sequences contain the motif\n\n";
	
	if ($negpresent > 0)
	{
		$count = 0;
		for (my $i=0; $i<@$negset; $i++)
		{
			my $target = ${$negset}[$i];
			my @positions = ();
			my @strings = ();
			if ($containsgap < 1)
			{
				while ($target =~ /$pattern/g)
				{
					my $nextpos = pos $target; # pos points to the next character
					$nextpos = $nextpos-$patternlength+1; # +1 to start counting at 1 instead of 0
					push(@positions,$nextpos);
					my $motif = substr($target,$nextpos,$patternlength);
					push(@strings,$motif);
					pos $target = $nextpos + 1;
				}
				if (@positions > 0)
				{
					my $seqid = $i+1; # to start counting at 1 instead of 0
					print OUT "${$negheaders}[$i] (sequence $seqid)\n";
					print OUT "    position(s): @positions\n";
					print OUT "    motif(s): @strings\n";
					$count++;
				}
			}
			else
			{
				my $found = 0;
				if ($target =~ /$pattern/)
				{
					$found=1;
				}
				if ($found > 0)
				{
					my $seqid = $i+1; # to start counting at 1 instead of 0
					print OUT "${$negheaders}[$i] (sequence $seqid)\n";
					$count++;
				}
			}
		}
		print OUT "$count negative sequences contain the motif\n\n\n";
	}
}

# check the coverage for the ensemble of motifs in the input file
sub check_coverage()
{
	open(IN,"$outputfile");
	while ($line=<IN>)
	{
		if ($line =~ />\s*(\S+)/)
		{
			my $name = $1;
			if (exists($coverage_hash{$name}))
			{
				$coverage_hash{$name} += 1;
			}
			else
			{
				$coverage_hash{$name} = 1;
			}
		}
	}
	close(IN);
	open(OUT,">>$outputfile");
	print OUT "\n\n************************************\n\n";
	print OUT "COVERAGE\n";
	print OUT "********\n";
	
	foreach my $value (sort {$coverage_hash{$b} <=> $coverage_hash{$a}} keys(%coverage_hash))
	{
		print OUT "$value ($coverage_hash{$value} motifs match)\n";
	}
	
	print OUT "\n" . keys(%coverage_hash) . " sequences covered in total\n";
}

# load a classification scheme defined by the user (see manual)
sub load_user_hierarchy()
{
	my $file = $_[0];
	open(IN,$file) or die "Can not open hierarchy file $file\n";
	my $line = <IN>;
	while ($line !~ /Definitions/)
	{
		$line = <IN>;
	}
	$line = <IN>;
	while ($line =~ /=/)
	{
		chomp($line);
		$line =~ /(\S+)\s*=\s*(.+)/ or die "Wrong line format in hierarchy file: $line\n";
		my $class = $1;
		my $values = $2;
		$values =~ s/ //g; #remove spaces
		$values =~ s/,//g; #remove commas
		$string_representation{$class} = "\[$values\]";
		$line = <IN>;
	}
}

# checks if pattern contains gap
sub contains_gap()
{
	my $patternstring = $_[0];
	if ($patternstring =~ /\sgap\s/)
	{
		return 1;
	}
	return 0;	
}
