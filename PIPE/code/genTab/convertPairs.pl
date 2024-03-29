#!/usr/bin/perl

$SEQ_FILE="data/protein_sequences.txt";
$STR_PAIR_FILE="data/protein_pairs.txt";
$IDX_PAIR_FILE="data/protein_pairs_index.txt";

open(SEQFILE, $SEQ_FILE);
$numProteins=0;
%nameHash=();
while(<SEQFILE>) {
	chomp;
	($name, $seq) = split(/\s+/);
	if (($name ne "") && ($nameHash{$name} eq "")) {
		$nameArray[$numProteins] = $name;
		$nameHash{$name} = $numProteins++;
	}
}
close(SEQFILE);


open(PAIRFILE, $STR_PAIR_FILE);
@idxArray=();
while(<PAIRFILE>) {
	chomp;
	($p1, $p2) = split(/\s+/);
	if ($nameHash{$p1} eq "") {
		print "ERROR: $p1 is not in sequence file\n";
		exit -1;
	}
	if ($nameHash{$p2} eq "") {
		print "ERROR: $p2 is not in sequence file\n";
		exit -1;
	}
	$p1Idx = $nameHash{$p1};
	$p2Idx = $nameHash{$p2};

	if ($idxArray[$p1Idx] eq "") {
		$idxArray[$p1Idx] = "$p2Idx";
	}
	else {
		$idxArray[$p1Idx] .= ":$p2Idx";
	}

	if ($idxArray[$p2Idx] eq "") {
		$idxArray[$p2Idx] = "$p1Idx";
	}
	else {
		$idxArray[$p2Idx] .= ":$p1Idx";
	}

}
close(PAIRFILE);


open(IDXPAIRFILE, ">$IDX_PAIR_FILE");
print IDXPAIRFILE "$numProteins\n";
$cnt=0;
for ($i=0; $i<$numProteins; $i++) {
	%temp = ();
	@list = grep ++$temp{$_} < 2, split(":", $idxArray[$i]);
	$numNeighbours = scalar(@list);
	print IDXPAIRFILE "$nameArray[$i] $numNeighbours @list\n";
	$cnt += $numNeighbours;
}
close(IDXPAIRFILE);

print "Average Neighbours per Protein: " . $cnt/$numProteins . "\n";
