#!/usr/bin/ruby

class Protein
	attr_reader :accession, :ncbi, :description, :sequence

	def initialize(accession, ncbi, description, sequence)
		@accession = accession
		@ncbi = ncbi
		@description = description
		@sequence = sequence
	end
end

class ProteinDb
	def initialize(filename)
		@filename = filename

 		# hash from acid sub-sequence (eg. "EQTW") to protein accession
		@index = {}

		# the length of the subsequences we index
		@index_length = 4

		# the protein db loaded in from the file: a hash from accession to 
		# Protein
		@db = {}

		load_database
	end

	# parse a header line from a fasta db
	# cases:
	# gi|136429|sp|P00761|TRYP_PIG TRYPSIN P
	# gi|4758488|ref|NP_004119.1| general tr
	# gi|229552|prf||754920A albumin [Bos pr
	def parse_header(header)
		if not header =~ /gi\|(\d+)\|(.*)\|(...\d+(.\d+)?)?\| ?(.*)/
			return nil
		end
	
		accession = $~[1].to_i
		ncbi = $~[3]
		description = $~[5]
	
		return accession, ncbi, description
	end

	# load a protein name database
	def load_database
		log "loading #{@filename} ..."

		n_errors = 0
		n_predicted = 0
		mismatches = []

		file = ReadFile.new(@filename)

		while file.line do
			accession, ncbi, description = parse_header(file.line)

			if accession == nil
				mismatches << file.line

				# skip to next ">" line
				while file.next and file.line[0] != 62 do
				end
				next
			end

			# skip predicted but unverified proteins
			if description =~ /^PREDICTED: /
				n_predicted += 1

				# skip to next ">" line
				while file.next and file.line[0] != 62 do
				end
				next
			end

			# we can drop any " [Homo sapiens]" at the end
			# there are some neanderthalensis ones in there, leave them
			if description =~ /(.*) \[Homo sapiens\]/
				description = $~[1]
			end

			# read the sequence
			sequence = ""
			while file.next and file.line[0] != 62 do
				sequence += file.line
			end

# 			puts file.line
# 			puts "\taccession = #{accession}"
# 			puts "\tncbi = #{ncbi}"
# 			puts "\tdescription = #{description}"
# 			puts "\tsequence = #{sequence}"

			if @db.has_key?(accession)
				err "accession collision! #{accession} appears twice"
				next
			end

			@db[accession] = Protein.new accession, ncbi, description, sequence
		end

		log "loaded #{@db.length} proteins from #{@filename}"

		if n_predicted > 0 
			log "ignored  #{n_predicted} 'PREDICTED:' proteins"
		end

		if mismatches != []
			log "#{mismatches.length} protein DB mismatches:"
			mismatches.each {|mismatch| log mismatch}
			log "(end of mismatch list)"
		end
	end

	# build the sub-sequence index
	def build_index
		if @index != {}
			return
		end

		log "indexing protein database ..."
	
		@db.each do |accession, protein|
			sequence = protein.sequence
			(0 ... sequence.length - @index_length - 1).each do |i|
				subseq = sequence[i .. i + @index_length - 1]
	
				if not @index.has_key?(subseq)
					@index[subseq] = []
				end
				@index[subseq] << accession
			end
		end
	
	#	@index.each do |key, value| 
	#		log "#{key} in #{value.length} proteins"
	#	end
	end

	# lookup a single protein by accession
	def by_accession(i)
		@db[i]
	end

	# find a list of all accession numbers of prots which have or contain the 
	# sequence
	def by_sequence(sequence)
		hits = []

		build_index

		# look up the first 4 acids in the index, just search those
		subseq = sequence[0 .. @index_length - 1]

		if not @index.has_key?(subseq)
			err "subsequence #{subseq} not found in index!"
		else
			@index[subseq].each do |accession|
				protein = by_accession(accession)
				if protein.sequence =~ /#{sequence}/
					hits << protein.accession
				end
			end
		end
	
		return hits
	end
end
