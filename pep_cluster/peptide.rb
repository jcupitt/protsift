#!/usr/bin/ruby

# a peptide from Odu's CSV files

class Peptide
	attr_reader :sample, :region, :sequence, :filename, 
		:charge, :score, :norm, :raw, :filename

	def initialize(sample, region, sequence, filename, 
		charge, score, norm, raw)
		@sample = sample
		@region = region
		@sequence = sequence
		@filename = filename
		@charge = charge
		@score = score
		@norm = norm
		@raw = raw
	end
end

class PeptideList
	include Enumerable

	# a threshold for each charge
	# xcorr_threshold = [2, 3.2, 3.4]
	@@xcorr_threshold = [1.5, 2.0, 2.5]

	# remove peps with deltacn below this
	@@deltacn_threshold = 0.2

	# remove prots with fewer than this many distinct peps
	@@pep_thresh = 2

	def initialize(db, filename, &filter_block)
		@db = db
		@filename = filename
		@filter_block = filter_block

		# this is the thing we build: a hash from accession to sequence to 
		# Peptide
		@table = Hash.new {|h, k| h[k] = Hash.new(&h.default_proc)}

		load_file
	end

	def load_file
		srf = Ms::Sequest::Srf.new(@filename)
		
		n_dupes = 0
		n_nils = 0
		n_deltacn = 0
		n_filter = 0
		n_xc = 0

		srf.out_files.each do |outfile|
			# just take the first hit from each file
			# we use deltacn to throw away peps from files with two very close
			# hits
			pep = outfile.hits.first

			# a few files in every SRF seem to have no hits, they have some kind
			# of metadata in there instead
			if pep == nil
				n_nils += 1
				next
			end

			if pep.charge < 1 or pep.charge > 3
				err "argh! pep.charge = #{pep.charge}"
				return 
			end

			# Xc too low
			if pep.xcorr < @@xcorr_threshold[pep.charge - 1]
				n_xc += 1
				next
			end

			# deltacn too low
			if pep.deltacn < @@deltacn_threshold
				n_deltacn += 1
				next
			end

			sequence = pep.aaseq

			# find all proteins containing this sequence
			all_prots = @db.by_sequence(sequence)

			if all_prots == []
				err "pep #{sequence} not found in database"
				next
			end

			# filter proteins with our block, if any
			if @filter_block
				all_prots.delete_if do |accession|
					protein = @db.by_accession(accession)

					test = @filter_block.call(protein)

					if test
						n_filter += 1
					end

					test
				end
			end

			# avoid making the pep if we can
			if all_prots != []
				peptide = Peptide.new sequence, pep.xcorr, pep.charge, pep.mh, 
					pep.first_scan, pep.last_scan, @filename
			end

			all_prots.each do |accession|
				# is the peptide there already? we want to filter out 
				# duplicates 
				# if it is there already, don't replace if the one that's 
				# there now has a better xc
				if @table[accession].has_key?(sequence) 
					n_dupes += 1

					if @table[accession][sequence].xc > pep.xcorr
						next
					end
				end

				@table[accession][sequence] = peptide
			end
		end

		if n_nils > 0
			log "#{n_nils} dta files have no peps"
		end

		log "#{srf.out_files.length - n_nils} peptides"

		if n_dupes > 0
			log "#{n_dupes} duplicate peptides found"
		end

		if n_xc > 0
			log "#{n_xc} peptides failed Xc filtering"
		end

		if n_deltacn > 0
			log "#{n_deltacn} peptides failed deltaCn filtering"
		end

		if n_filter > 0
			log "#{n_filter} proteins failed filtering"
		end

		log "#{@table.length} proteins found"

		# now remove proteins with less than x peptides
		@table.delete_if {|accession, peptides| peptides.length < @@pep_thresh}

		log "#{@table.length} proteins after #{@@pep_thresh} pep filtering"
	end

	def each(&block)
		@table.each do |key, value|
			yield key, value
		end
	end
end
