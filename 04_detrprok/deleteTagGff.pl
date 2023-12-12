#!/usr/bin/perl -w

use vars qw($USAGE);
use strict;

=head1 NAME

deleteTagGff.pl - delete command-line specified tags of a gff file

=head1 SYNOPSIS

% deleteTagGff.pl -i gffFile [-d tag] [-h] 

=head1 DESCRIPTION
This script will parse artemis gff file into smart gff file
   -i|--input              fileName  input file name
   [-d|--delete tagId]     tag name (or list of tag names comma-separated) to delete
   [-h|--help]             help mode then die 

if "ID" or "Name" tags don't exist, create a minimal "Name" tag with value Name=column1+/-begin..end

=head1 AUTHOR - Claire Toffano-Nioche - apr.11

contexte d'usage : 
blast result analyses of Vibrio ncRNA candidates

=cut
my ($inFile, @delTags);

foreach my $num (0 .. $#ARGV) {
	SWITCH: for ($ARGV[$num]) {
	/--input|-i/	 && do { $inFile=$ARGV[$num+1] ; open (INFILE,"<$ARGV[$num+1]") or die "Can't open \"$ARGV[$num+1]\"\n" ;last; };
	/--delete|-d/    && do { @delTags=(@delTags,split(',',$ARGV[$num+1])); last; }; # push for multiple -d, split for comma-separated list of tags
	/--help|-h/      && do { exec("pod2text $0\n") ; die };
	}
}
while (my $line = <INFILE>) {
    chomp($line) ;
	my @info = split("\t", $line) ; 
	my @tags = split(";", $info[8]);
	#print join("|", @delTags), "\n" ;
	#print join("|", @tags), "\n" ;
	#print ":::::$#tags tags\n";
	my $defaultName="Name=".$info[0].$info[6].$info[3]."..".$info[4] ;
	# tags suppression:
	if ($delTags[0] ne "") { # if not null to delete tags only if -d used
		for (my $t = 0 ; $t <= $#delTags ; $t++) {
		    	#print "\ntag to delete : $delTags[$t]\n" ;
			my $i = 0;
			#my $supprOk = 0 ;
			while ($i <= $#tags) {
				if ($tags[$i] =~ /$delTags[$t]=/) {	
					@tags = @tags[0..($i-1),($i+1)..$#tags] ;  
			#		$supprOk = 1 ;
				#print join("|", @tags), "\n" ;
				#print "suppr $delTags[$t] ok:$#tags tags\n";
				}
				$i++
			}
		}
	}
	# tag Name treatment:
	my $IDexist = 0;
	my $i=0 ;
	#print "$#tags tags\n";
	while (($i <= $#tags) and (!$IDexist)) { 
		if (($tags[$i] =~ /Name/) 
		# or ($tags[$i] =~ /ID/)
		                         ) { # check for tags Name or ID
			$IDexist = 1 ;
			#print "exist $tags[$i]\n";
		} 
		$i++;
	}
	if (!$IDexist) {
		unshift(@tags, ($defaultName)) ; # add tag "Name=" if no Name tag
		#print "new tag name\n";
	}
	# result by display:
	my $col9=join(";", @tags);
	print join("\t", ($info[0], $info[1], $info[2], $info[3], $info[4], $info[5], $info[6], $info[7], $col9)), "\n" ;
}
close(INFILE) ;
exit(0) ;
