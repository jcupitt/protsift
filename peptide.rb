#!/usr/bin/ruby

# walk an srf file and return a hash from good-protein-reference to a hash
# from peptide-name to peptide
#
# a pep looks like this:
# 	<Ms::Sequest::Srf::Out::Pep 
# 		aaseq="RALRNNTDPYK" 
# 		sequence="R.RALRNNTDPYK.R"
# 		mh=1347.712772606 
# 		deltacn_orig=0.0 
# 		sf=0.000485505705000833
# 		sp=1.03858149051666 
# 		xcorr=0.479519844055176 
# 		id=18669 
# 		rsp=10 
# 		ions_matched=2
# 		ions_total=20 
# 		prots(#)=1 
# 		deltamass=-1.71398520650018 
# 		ppm=1270.15801085689
# 		base_name="P7R7_090402_10" 
# 		first_scan=764 
# 		last_scan=764 
# 		charge=1
# 		deltacn=0.00902968738228083 
# 		srf(base_name)="P7R7_090402_10" 
# 	>

# filter parameters

# we can't use the SRF pep class, it drags in too much stuff and we run out 
# of memory
class Peptide
	attr_reader :sequence, :sf, :xc, :charge, :mh, :first_scan, 
		:last_scan, :filename

	def initialize(sequence, sf, xc, charge, mh, 
				   first_scan, last_scan, filename)
		@sequence = sequence
		@sf = sf
		@xc = xc
		@charge = charge
		@mh = mh
		@first_scan = first_scan
		@last_scan = last_scan
		@filename = filename
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
		@broken = false

		# this is the thing we build: a hash from accession to sequence to 
		# Peptide
		@table = Hash.new {|h, k| h[k] = Hash.new(&h.default_proc)}

		load_file
	end

	def is_broken?
		@broken
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
				err "argh! pep.charge = #{pep.charge}, filename = #{@filename}"
				@broken = true

				# this file has been corrupted, bail out
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

			# filter matching proteins with our block, if any
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

			if all_prots == []
				# log "pep #{sequence} filtered"
				next
			end

			if pep.sf > 1.0 or pep.sf < 0.0001 
				log "pep = #{pep.inspect}"
			end

			peptide = Peptide.new sequence, pep.sf, pep.xcorr, pep.charge, 
					pep.mh, pep.first_scan, pep.last_scan, @filename

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

	def each
		@table.each do |key, value|
			yield key, value
		end
	end
end
