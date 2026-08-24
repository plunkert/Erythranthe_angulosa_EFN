Readme for "No evidence of extrafloral nectaries in Erythranthe angulosa" Martin-Eberhardt, Smith, & Plunkert. 

Data were collected from seed-grown *Erythranthe angulosa* plants growing in growth chambers at Michigan State University. 
Seeds were from a single ripe fruit ("accession 1") collected at the type locality for this species 
(36.467044 -118.134064, Nesom and Berger 2020) in early July 2024 and germinated from soil collected from beneath this same mother plant. 

Induction experiment: plants were germinated after 14 days of stratification at 4 C. Plants were grown in 14 / 10 light dark, 
28.3 C, 10% humidity during day, 15.5 C, 35% humidity night. Treatments (clipping, JA, or soil amendment) 
were imposed when plants were six weeks old and repeated every 2 weeks for a total of three months from first treatment to 
data collection. At the end of the induction experiment, we collected data on both the number of axillary swellings (hypothesize EFNs)
and morphological data (leaf size, stem number, stem thickness) on two distinct morphs germinated from soil and included in the 
experiment. All data were collected on Oct 2nd, 2025. 


Data:
	Induction_exp.csv -- data from an EFN induction experiment applying three treatments (clipping, Jasmonic acid application, 
and silica soil amendment) and testing for a change in axillary swelling (hypothesized EFN) production
				or induction of nectar production.
	Phenotype_morphology.csv -- data collected on the morphology (stem thickness, leaf size, and stem number) of progeny from 
two distinct morphs observed germinating from the soil

Code:
	analyses.r -- all data cleaning and analyses reported in the paper


Induction_exp.csv
	plant_num -- unique ID for that plant
	family -- the accession (genetic family) of the individual
			1-1G2 -- 'accession 1', seed from single ripe fruit of plant with axillary swellings
			11G, 12G, 50G, 51G -- families of the 'Robust Morph' -- numbers correspond to mother size, with 11 and 12 
					          being the smallest within this morph and 50 and 51 the largest
                        53G, 54G, 77G, 78G -- families of the 'Thin-stemmed Morph' -- numbers correspond to mother size
	treatment -- the treatment applied to the individual
			clip -- 1/3 of leaf clipped 
			JA -- 1mM jasmonic acid sprayed on leaves
			soil -- plants grown in a mixture of potting mix, crushed granite, and silica sand (See paper for details)	
			control -- plants in standard SureMix potting mix, not clipped or exposed to JA
	axillary_swellings -- number of axillary swellings observed on the entire plant
	notes -- notes

Phenotype_morphology.csv
	family -- see above
	treatment -- see above
	stem_thickness -- the diameter of the central stem (mm)
	len_lengest_leaf -- the length of the longest leaf (mm)
	stem_num -- the number of stems the plant had
	plant_num -- unique ID for that plant
	notes -- notes
			