#!/usr/bin/ruby

require 'csv'
require 'set'
require 'ftools'
require 'optparse'

# get the files from the same directory as us
$: << File.dirname(__FILE__) 

require 'read_file.rb'
require 'protein.rb'
require 'peptide.rb'

$options = {
	:verbose => true,
	:fatal => false,
	:output => "output.csv",
	:output2 => "output2.csv"
}
OptionParser.new do |opts|
	opts.banner = "Usage: #{$0} [options] row17.csv ..."

	opts.on("-o", "--output=FILE", "output to FILE") do |file|
		$options[:output] = file
	end
	opts.on("-f", "--fatal", "Stop on first error") do |v|
		$options[:fatal] = v
	end
	opts.on("-q", "--quiet", "Run quietly") do |v|
		$options[:verbose] = false
	end
end.parse!

def log(msg)
	if $options[:verbose]
		puts msg
	end
end

def err(msg)
	puts msg

	if $options[:fatal]
		exit
	end
end

# auto-expanding hash
peptide_table = Hash.new {|h, k| h[k] = Hash.new(&h.default_proc)}

n_peps = 0

ARGV.each do |filename|
	if filename !~ /row(.*)\.csv/
		err "bad filename #{filename}"
		next
	end
	sample = "xx"
	region = $~[1].to_i

	log "loading #{filename} ..."

	linen = 0
	CSV.open(filename, 'r') do |line|
		# all useful ones should begin with a "g"
		if line[130] !~ /^g.*$/
			next
		end

		# column EA is index 130 and should be something like
		# "gi|61743954|ref|NP_001611.1|"
		if line[130] !~ /^gi\|(\d+)\|ref\|(.._\d+(\.\d+)?)\|(REVERSED)?$/
			err "mismatch for #{filename}, #{line[74]}"
			next
		end
		accession = $~[1].to_i
		ncbi = $~[2]
		reversed = $~[4]

		if reversed
			next
		end

		score = line[129].to_f
		sequence = line[131]
		old_description = line[133]

		charge = line[4].to_i

		norm = line[8..37]
		raw = line[38..67]
		intens = line[68..97]
		reten = line[98..127]

		# replace if we have a better score
		old_pep = peptide_table[sample][region][sequence]
		if old_pep == {} or old_pep.score < score
			peptide_table[sample][region][sequence] = 
				Peptide.new(sample, region, sequence, filename, 
						charge, score, norm, raw)
				n_peps += 1
		end
	end
end

log "found #{n_peps} peptides"

# load the human protein db
protein_db = ProteinDb.new "human_protein2010.fasta"

log "building protein table ..."

# build a table of proteins, with a set of peps for each one
protein_table = Hash.new {|h, k| h[k] = Set.new}

peptide_table.each do |sample, regions|
	regions.each do |region, sequences|
		sequences.each do |sequence, peptide|
			protein_db.by_sequence(peptide.sequence).each do |accession|
				protein_table[accession] << peptide
			end
		end
	end
end

log "found #{protein_table.length} proteins"

log "selecting and combining regions for each protein ..."

protein_table.each do |accession, peptides|
	hits = Hash.new(0)
	peptides.each do |peptide|
		hits[peptide.region] += 1
	end
	best_region, n_hits = hits.max
	peptides.delete_if {|peptide| (peptide.region - best_region).abs > 1}
end

log "removing duplicate peptides ..."

n_peps = 0
protein_table.each {|accession, peptides| n_peps += peptides.length}
log "(#{n_peps} peptides before dedupe)"

protein_table.each do |accession, peptides|
	unique = {}
    peptides.each do |peptide|
		old_pep = unique[peptide.sequence]
		if old_pep == nil or
			old_pep.score < peptide.score
			unique[peptide.sequence] = peptide
		end
	end

	pep_array = Set.new
	unique.each {|sequence, peptide| pep_array << peptide}
	protein_table[accession] = pep_array
end

n_peps = 0
protein_table.each {|accession, peptides| n_peps += peptides.length}
log "(#{n_peps} peptides after dedupe)"

log "removing proteins with less than two peps ..."

protein_table.delete_if {|accession, peptides| peptides.length < 2}

log "(#{protein_table.length} proteins remain)"

log "building protein groups ..."

# we pair an array of accession numbers with a set of peptides ... those
# proteins could all be identified by that set of peps
protein_group = []

protein_table.each do |accession, peptides|
	# does a protein previously added to protein_group already have either a
	# subset or superset of this set of peptides?
	found = false

	protein_group.each_index do |i|
		accessions, peptides_group = protein_group[i]

		if peptides_group.subset?(peptides) or
			peptides.subset?(peptides_group) 
			# yes, add to this protein group
			protein_group[i][0] << accession

			# if we have a superset of the peps already recorded, swap our
			# larger set of peps for the one there already
			if peptides.length > peptides_group.length
				protein_group[i][1] = peptides
			end

			found = true
			break
		end
	end

	# no existing group is either a super or subset, make a new group
	if not found
		protein_group << [[accession], peptides]
	end
end

log "(#{protein_group.length} protein groups found)"

log "linking peptides to groups ..."

# hash from peptide sequence to protein group
peptide_index = Hash.new {|h, k| h[k] = []}

protein_group.each do |group|
	group_accession, group_peptides = group

	group_peptides.each do |peptide|
		peptide_index[peptide.sequence] << group
	end
end

log "finding cross-group peptides ..."
peptide_index.delete_if {|sequence, groups| groups.length < 2}

log "(#{peptide_index.length} cross-group peptides found)"

log "linking groups ..."

clusters = []

peptide_index.each do |sequence, groups|
	# remove all clusters we hit
	hits = []
	groups.each do |group|
		clusters.each do |cluster| 
			if cluster.member?(group)
				hits << cluster 
			end
		end
	end
	clusters -= hits

	# join hits into a single large cluster
	join = []
	hits.each {|x| join += x}

	# and make sure all the groups that sequence hits are in there
	join |= groups

	clusters << join
end

log "(#{clusters.length} protein clusters found)"

# find all groups not in any cluster and add them as singletons
log "finding all unused groups ..."

unused = []

protein_group.each do |group|
	found = false

	clusters.each do |cluster|
		if cluster.member?(group)
			found = true
			break
		end
	end

	if not found
		unused << [group]
	end
end

log "(#{unused.length} unused groups)"

clusters += unused

log "writing #{$options[:output]} ..."
CSV.open($options[:output], 'w') do |write|
	row = []
	row << "accession"
	row << "ncbi"
	row << "description"
	write << row

	row = []
	row << ""
	row << "sequence"
	row << "charge"
	row << "score"
	row << "sample"
	row << "region"
	row << "norm"
	29.times {row << ""}
	row << "raw"
	29.times {row << ""}
	write << row

	clusters.each do |cluster|
		if cluster.length > 1
			row = []
			row << "** cluster of #{cluster.length} proteins"
			write << row
		end

		cluster.each do |group|
			accessions, peptides = group

			accessions.each do |accession|
				protein = protein_db.by_accession(accession)
				row = []
				row << protein.accession
				row << protein.ncbi
				row << protein.description
				write << row
			end
	
			row = []
			6.times {row << ""}
	
			# sum the norm and raw peps
			norm = Array.new(30){0}
			raw = Array.new(30){0}
			peptides.each do |peptide|
				peptide.norm.each_index {|i| norm[i] += peptide.norm[i].to_f}
				peptide.raw.each_index {|i| raw[i] += peptide.raw[i].to_f}
			end
	
			norm.each {|x| row << x}
			raw.each {|x| row << x}
	
			write << row
	
			peptides.each do |peptide|
				row = []
				row << ""
				row << peptide.sequence
				row << peptide.charge
				row << peptide.score
				row << peptide.sample
				row << peptide.region
				peptide.norm.each {|x| row << x}
				peptide.raw.each {|x| row << x}
				write << row
			end
		end
	end
end

log "writing #{$options[:output2]} ..."
CSV.open($options[:output2], 'w') do |write|
	row = []
	row << "cluster size"
	row << "group size"
	row << "accession"
	row << "ncbi"
	row << "description"
	3.times {row << ""}
	row << "norm"
	15.times {row << ""}
	row << "raw"
	15.times {row << ""}
	write << row

	clusters.each do |cluster|
		# find the cluster with the largest number of peps
		max_peps = 0
		max_group = nil
		cluster.each do |group|
			accessions, peptides = group

			if peptides.length > max_peps
				max_peps = peptides.length
				max_group = group
			end
		end
		accessions, peptides = max_group
		protein = protein_db.by_accession(accessions.first)

		# sum the norm and raw peps
		norm = Array.new(30){0}
		raw = Array.new(30){0}
		peptides.each do |peptide|
			peptide.norm.each_index {|i| norm[i] += peptide.norm[i].to_f}
			peptide.raw.each_index {|i| raw[i] += peptide.raw[i].to_f}
		end

		row = []
		row << cluster.length
		row << accessions.length
		row << protein.accession
		row << protein.ncbi
		row << protein.description

		3.times {row << ""}

		norm.each {|x| row << x}
		raw.each {|x| row << x}

		write << row
	end
end

