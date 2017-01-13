#!/usr/bin/perl

$keep_lines = @ARGV[0];

open my $handle, '<', $keep_lines;
chomp(my @keep = <$handle>);
close $handle;

$n = 1;
$k = 0;


while(<STDIN>) {

	$line = $_;
	
	if( $n == @keep[$k] ) {
		print $line;
		$k++;
	}

	$n++;
}