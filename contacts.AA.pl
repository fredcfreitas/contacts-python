#!/usr/bin/perl -w
################################################################################
# Contacts.pl is a simple script written to analyze All-Atom structure-based   #
# simulations run in Gromacs.   You must supply a contact file, with 6-12      #
# parameters, the .tpr file used for the simulation, and a trajectory file.    #
# You must also have Gromacs 4.0 installed on your machine                     #
# This is a very straightforward script that you can modify in any way you see #
# fit.  This software comes with absolutely no guarantee, but let us know if   #
# believe you found a bug.                                                     #
# Written by Paul Whitford, 10/26/08                                           #
# Minor modification  1/12/15                                           #
################################################################################

# Edit the settings so they fit your needs...

# GMXPATH is the path to your Gromacs executables
$GMXPATH="/home/pwhitfor/BIN/gromacs-4.0/bin/";


# CONTFILE is the file that defines the contacts.  Specific formatting must be
# followed: Copy the "pairs" terms from your AA Structure-based topology file.
# remove the "1" and reformat each line so it is space delimited.  If you have
# a blank line in the contact file, the program will probably crash.
# example formatting can be found at http://smog-server.org/contact_ex
$CONTFILE="somefilename";

# This program will calculate 2 versions of Q.  First, it will calculate all of
# the residue pairs that have at least 1 native contact formed (AAnew).  It will
# also calculate the number of atom-atom contacts formed (AA).
# CUTOFF is how you determine a contact is formed, or broken.  CUTOFF=1.5 means
# a native contact is considered formed if the involved atoms are within 150%
# of their native distance.
$CUTOFF=1.5;

# If you would only like to calculate the Q of every Nth frame, then change SKIPFRAMES to the desired frequency of output.

$SKIPFRAMES=1;

# End of settings



# this reads in contact list

open(CONaa,"$CONTFILE")  or die "no AA contacts file";
$CONNUMaa=0;
while(<CONaa>){
	$ConI=$_;
	chomp($ConI);
	$CONNUMaa ++;

	@TMP=split(" ",$ConI);
	$Iaa[$CONNUMaa]=$TMP[0];
	$Jaa[$CONNUMaa]=$TMP[1];
	# this gets the native distances from the 6-12 parameters and converts to Ang.
	$Raa[$CONNUMaa]=(2*$TMP[3]/$TMP[2])**(1/6.0)*10.0;
}
close(CONaa);


print "What is the name of the tpr file?\n";
$TPR=<STDIN>;
chomp($TPR);
# go through all of the xtc files
print "What  xtc (or trr) file would you like to analyze?\n";
$XTCfile=<STDIN>;
chomp($XTCfile);
print "\n\n\n\n\n";
print "This script will calculate the number of native contacts for frames n*$SKIPFRAMES (n starts at 0) for trajectory file $XTCfile.\n";
print "IMPORTANT: A contact is defined as any native pair listed in $CONTFILE that is withing $CUTOFF times the distance given in $CONTFILE\n";
print "\n\n\n\n\n";

# you can change the flags sent to gromacs if you want to reduce the Q sampling
# i.e. every other frame is analyzed.  The trajectory file is converted to a pdb
# and then this script analyzes the pdb.
`echo 0 | $GMXPATH/trjconv  -skip $SKIPFRAMES -s $TPR -o tmp.pdb -f  $XTCfile`;

# this will give the number of newQ contacts at each snapshot
open(AAnQ,">$XTCfile.AA.nQ") or die "couldn't open nQ file";

# this will give a list at each frame of which contacts are formed
# the numbering will correspond to the numbering in $XTCfile.Qn.index
open(AAnQI,">$XTCfile.AA.nQi") or die "couldn't open nQi file";

# these are for the atom-atom contacts
open(AAQ,">$XTCfile.AA.Q") or die "couldn't open Q file";
open(AAQI,">$XTCfile.AA.Qi") or die "couldn't open Qi file";



open(PDB,"tmp.pdb") or die "pdb file missing somehow";
#GRAB THE ATOM residue numbers
while(<PDB>){
     	$LINE=$_;
       	chomp($LINE);

        if(substr($LINE,0,5) eq "MODEL"){
        	# start the residue storing
               	$ATOMNUM=0;
        }

        if(substr($LINE,0,4) eq "ATOM"){
           	$ATOMNUM=$ATOMNUM+1;
		$RES[$ATOMNUM]=substr($LINE,22,4);
        }

        if(substr($LINE,0,3) eq "TER"){
		last;
	}

}

close PDB;



# find the unique residue contacts and index them
$NEWcn=0;

for($i=1;$i<=$CONNUMaa;$i++){
	#see if the residue pair is a contact already.
	$UNIQUE=1;
	for($j=1;$j<=$NEWcn;$j++){
	# if it already appears, then change $UNIQUE
		if(($NEWCONi[$j] == $RES[$Iaa[$i]]) && ($NEWCONj[$j] ==$RES[$Jaa[$i]])){
			$UNIQUE=0;
			$MATCH=$j;
		}
	}

	if($UNIQUE){
		# if it is unique, then add it to the list
		$NEWcn++;
		$NEWCONi[$NEWcn]=$RES[$Iaa[$i]];
       		$NEWCONj[$NEWcn]=$RES[$Jaa[$i]];
       		$ATOMtoRESContNUM[$i]=$NEWcn;
		$NEWcount[$NEWcn]=1;
	}else{
		$ATOMtoRESContNUM[$i]=$MATCH;
		$NEWcount[$MATCH]++;
	}
}


open(OUT1,">$XTCfile.nQ.index") or die "something wrong happened";
print OUT1 "This the contact list refered to in $XTCfile.AA.nQi (residue i, residue j of contact 1)\n";
for($i=1;$i<=$NEWcn;$i++){
	print OUT1 "$NEWCONi[$i] $NEWCONj[$i]\n"
}

print OUT1 "This is the mapping from AA contacts to residue contacts\n";
for($i=1;$i<=$CONNUMaa;$i++){
        print OUT1 "$ATOMtoRESContNUM[$i] \n"
}
close OUT1;





open(PDB,"tmp.pdb") or die "pdb file missing somehow";
###$SAMPLES=0; #I think this is obsolete
	#GRAB THE ATOM COORDS
while(<PDB>){
	$LINE=$_;
	chomp($LINE);

	if(substr($LINE,0,5) eq "MODEL"){
		$ATOMNUM=0;
	}


        if(substr($LINE,0,4) eq "ATOM"){
	        # store positions, index and residue number
		$ATOMNUM=$ATOMNUM+1;
	        $X[$ATOMNUM]=substr($LINE,30,8);
                $Y[$ATOMNUM]=substr($LINE,38,8);
               	$Z[$ATOMNUM]=substr($LINE,46,8);
        }


        if(substr($LINE,0,3) eq "TER"){

	# do contact calculations

	# determine the newQ contacts
               	for($i=1;$i<=$NEWcn;$i++){
                       	$FORMEDnQ[$i]=0;
                }

                for($i = 1;$i <= $CONNUMaa; $i ++){
                        $R=sqrt( ($X[$Iaa[$i]]-$X[$Jaa[$i]])**2.0 +
                                 ($Y[$Iaa[$i]]-$Y[$Jaa[$i]])**2.0 +
                                 ($Z[$Iaa[$i]]-$Z[$Jaa[$i]])**2.0);

                       	if($R < $CUTOFF*$Raa[$i]){
				$FORMEDnQ[$ATOMtoRESContNUM[$i]]=1;
                        }
                }

		$Qaa=0;
		for($i=1;$i<=$NEWcn;$i++){
			if($FORMEDnQ[$i]){
				$Qaa++;
                                print AAnQI "$i\n";
			}

		}
                print AAnQ "$Qaa\n";

		# determine the atom-atom contacts

               	$Qaa=0;

		for($i = 1;$i <= $CONNUMaa; $i +=1){
			$R=sqrt( ($X[$Iaa[$i]]-$X[$Jaa[$i]])**2.0 +
			       	 ($Y[$Iaa[$i]]-$Y[$Jaa[$i]])**2.0 +
  			 	 ($Z[$Iaa[$i]]-$Z[$Jaa[$i]])**2.0);

 			if($R < $CUTOFF*$Raa[$i]){
				$Qaa ++;
				print AAQI "$i\n";
			}
		}
	        print AAQ "$Qaa\n";
        }
}
close(AAQ);
close(AAQI);

`rm tmp.pdb`;
