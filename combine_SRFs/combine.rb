#!/usr/bin/ruby

require 'csv'
require 'set'
require 'ftools'
require 'optparse'

require 'rubygems'
require 'ms/sequest/srf'

# get the files from the same directory as us
$: << File.dirname(__FILE__) 

require 'read_file.rb'
require 'protein.rb'
require 'fraction.rb'
require 'xme_filter.rb'
require 'rename.rb'
require 'peptide.rb'

$options = {
	:verbose => true,
	:fatal => false,
	:output => "output.csv",
	:output2 => "output2.csv"
}
OptionParser.new do |opts|
	opts.banner = "Usage: #{$0} [options] NP4Row10feature.csv ..."

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

rename_db = RenameDb.new "protein_rename"
protein_db = ProteinDb.new "human_protein2010.fasta"
xme_list = XMEFilter.new protein_db, rename_db, "All_XMEs_Found.csv"
# xme_list.save_database($options[:output])
fraction_db = Fraction.new "protein_fraction.csv"

# now trim the protein DB down to just the XMEs
n_xmes = 0
protein_db.delete_if do |protein|
	if xme_list.find_desc(protein.description)
		n_xmes += 1
		false
	else
		true
	end
end

log "(#{n_xmes} prots survive xme filtering)"

# places we search for srfs
srfs = [
	"3DSC+MTT/Cyto/All",
	"3DSC+MTT/Micro/All",
	"FemaleVSLiver/Exp2_micro/All",
	"FemaleVSLiver/Exp3_cyto/All"

	#"MTT_micro/All",
	#"MTT_cyto/All"
]

# sample names we know
sample_names = [
	"1188",
	"254",
	"EPISKIN",
	"RHE",
	"HC",
	"D",
	"L",
]

# build a protein table here
# a big nested hash 
# protein_table[desc][sample] = list of peptides
protein_table = Hash.new {|h, k| h[k] = Hash.new(&h.default_proc)}

n_files = 0
srfs.each do |base|
	Dir.glob("/home/sven/" + base + "/*.srf").each do |filename|
		if File.basename(filename, '.srf') =~ /(.*)[-_]([ABCDEFGHIJKLMNOPQRSUVWXYZ0-9]*)(T\d+)?([a-z])?/
			n_files += 1
		end
	end
end

log "#{n_files} .srf files to scan"
log "estimate #{(n_files / 60).to_f} minutes to process"

srfs.each do |base|
	log "scanning #{base} ..."
	Dir.glob("/home/sven/" + base + "/*.srf").each do |filename|
		if not File.basename(filename, '.srf') =~ /(.*)[-_]([ABCDEFGHIJKLMNOPQRSUVWXYZ0-9]*)(T\d+)?([a-z])?/
			err "bad filename #{filename}\n"
			next
		end
		region = $~[1].to_i
		sample = $~[2]
		timepoint = $~[3]
		replicate = $~[4]
		number = nil

		if filename =~ /yto/
			fraction = "C"
		else
			fraction = "M"
		end

		if sample =~ /D(\d+)/
			sample = "D"
			number = $~[1].to_i
		end

		if sample =~ /L(\d+)/
			sample = "L"
			number = $~[1].to_i
		end

		if not sample_names.include?(sample)
			err "uknown sample #{sample}\n"
			next
		end

		# log "filename = #{filename}"
		# log "region = #{region}"
		# log "sample = #{sample}"
		# log "timepoint = #{timepoint}"
		# log "replicate = #{replicate}"
		# log "number = #{number}"

		log "processing #{sample}, region #{region} ..."

		peptide_list = PeptideList.new protein_db, filename 

		if peptide_list.is_broken?
			log "broken: #{filename}"
			next
		end

		# merge with our global table
		peptide_list.each do |accession, peptides|
			peptides.each do |sequence, peptide|
				protein = protein_db.by_accession(accession)
				desc = xme_list.find_desc(protein.description)
				sequence = fraction + "-" + peptide.sequence
				xc = peptide.xc

				if not protein_table[desc][sequence].has_key?(sample) or
					protein_table[desc][sequence][sample].xc < xc
					protein_table[desc][sequence][sample] = peptide
				end
			end
		end
	end
end

filename = $options[:output]
log "writing #{filename} ..."
CSV.open(filename, 'w') do |write|
	row = []
	row << ""
	row << ""
	row << ""
	sample_names.each {|sample_name| row << sample_name}
	write << row

	row = []
	row << "Protein"
	write << row

	row = []
	row << ""
	row << "Fraction"
	row << ""
	row << "Sample max Sf"
	write << row

	row = []
	row << ""
	row << ""
	row << ""
	row << "Sample npeps"
	write << row

	row = []
	row << ""
	row << "Sequence"
	row << "Charge"
	row << "Sample max Xc"
	6.times {row << ""}
	row << "Peptide also appears in"
	write << row

	protein_table.keys.sort.each do |desc|
		sequences = protein_table[desc]
		fraction = fraction_db[desc]

		found = false
		sequences.each do |sequence, samples|
			if sequence =~ /#{fraction}-[A-Z]+/
				found = true
				break
			end
		end
		if not found
			next
		end

		num_peps = Hash.new(0)
		max_sf = Hash.new(0)
		charges = Hash.new(0)
		sequences.each do |sequence, samples|
			if sequence =~ /#{fraction}-[A-Z]+/
				sample_names.each do |sample_name|
					if samples.has_key?(sample_name)
						num_peps[sample_name] += 1
						max_sf[sample_name] = [max_sf[sample_name], samples[sample_name].sf].max
						charges[sequence] = samples[sample_name].charge
					end
				end
			end
		end

		row = []
		row << desc + ", " + xme_list[desc].to_a.join(", ")
		write << row

		row = []
		row << ""
		row << fraction
		row << ""
		sample_names.each {|x| row << "%.3f" % max_sf[x]}
		write << row

		sequences.keys.sort.each do |sequence|
			samples = sequences[sequence]

			sequence =~ /([CM])-([A-Z]+)/
			frac = $~[1]
			seq = $~[2]

			if frac != fraction
				next
			end

			row = []
			row << ""
			row << seq
			row << charges[sequence]

			sample_names.each do |sample_name|
				if samples.has_key?(sample_name)
					row << "%.2f" % samples[sample_name].xc
				else
					row << "-"
				end
			end

			hits = protein_db.by_sequence(seq)
			prots = hits.map {|x| protein_db.by_accession(x)}
			prots = prots.map {|x| x.ncbi}
			prots.delete_if {|x| xme_list[desc].include?(x)}
			if prots.length > 0 
				row << prots.join(", ")
			else
				row << "-"
			end

			write << row
		end

		row = []
		row << ""
		row << ""
		row << ""
		sample_names.each {|x| row << num_peps[x]}
		write << row
	end
end
