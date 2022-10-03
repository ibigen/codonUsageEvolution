'''
Created on Oct 31, 2017

@author: mmp
'''

import os, random, ntpath, sys, getpass, math, stat, gzip
from time import strftime, localtime
from Bio import SeqIO
from CAI import CAI
import Bio.Data.CodonTable as ct


class Utils(object):
	'''
	class docs
	'''
	TEMP_DIRECTORY = os.getenv("TMP", "/tmp")  ## can be of the user
	TEMP_DIRECTORY_OS = "/tmp"  ## it is always the OS temp diretory, in Cluster is SSD
	COUNT_DNA_TEMP_DIRECTORY = "extract_dna"
	ALL_aminos = "AllAminos"
	DICT_bases = {'A': 0, 'C': 1, 'T': 2, 'G': 3}

	def __init__(self, genetic_code=1):
		'''
		Constructor
		'''
		self.genetic_code = genetic_code

	def get_unique_file(self, file_name):
		"""
		get unique file name from a file_name
		return '<path file_name>/<random number>_<file_name>'
		"""
		temp_file_name = "{}_{}".format(random.randrange(10000000, 99999999, 10), ntpath.basename(file_name))
		main_path = os.path.dirname(file_name)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			if (not os.path.exists(os.path.join(main_path, temp_file_name))): break
			temp_file_name = "{}_{}".format(random.randrange(10000000, 99999999, 10), ntpath.basename(file_name))
		return os.path.join(main_path, temp_file_name.replace(" ", "_"))

	def get_temp_file(self, file_name, sz_type):
		"""
		return a temp file name
		"""
		main_path = os.path.join(self.TEMP_DIRECTORY, getpass.getuser(), self.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)):
			os.makedirs(main_path)
		else:
			temp_path = os.path.join(self.TEMP_DIRECTORY, getpass.getuser())
			if os.name != 'nt':
				cmd = "touch {}".format(temp_path)
				os.system(cmd)

		while 1:
			return_file = os.path.join(main_path, "cai_" + file_name + "_" + str(
				random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass

	def get_temp_file_from_dir(self, dir_out, file_name, sz_type):
		"""
		return a temp file name
		"""
		if (not os.path.exists(dir_out)): os.makedirs(dir_out)
		while 1:
			return_file = os.path.join(dir_out, "cai_" + file_name + "_" + str(
				random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass

	def get_temp_dir(self, b_os_temp_dir=False):
		"""
		:param b_os_temp_dir True if "/tmp" dir is mandatory, best approach for quick calls to file system (Cluster)
		return a temp directory
		"""
		main_path = os.path.join(self.TEMP_DIRECTORY_OS if b_os_temp_dir \
									 else self.TEMP_DIRECTORY, getpass.getuser(), self.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)):
			os.makedirs(main_path)
		else:
			temp_path = os.path.join(self.TEMP_DIRECTORY, getpass.getuser())
			if os.name != 'nt':
				cmd = "touch {}".format(temp_path)
				os.system(cmd)

		while 1:
			return_path = os.path.join(main_path, "cai_{}_{}".format(
				str(os.getpid()), str(random.randrange(10000000, 99999999, 10))))
			if (not os.path.exists(return_path)):
				os.makedirs(return_path)
				return return_path

	def get_file_name_without_extension(self, file_name):
		"""
		return file name without extension
		"""
		return os.path.splitext(os.path.basename(file_name))[0]

	def remove_temp_file(self, sz_file_name):
		"""
		prevent to remove files outside of temp directory
		"""
		if (sz_file_name == None): return

		if os.path.exists(sz_file_name) and len(sz_file_name) > 0 and sz_file_name.startswith(self.TEMP_DIRECTORY):
			cmd = "rm " + sz_file_name
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to remove a file")

	def remove_file(self, sz_file_name):
		"""
		Remove files
		return True if the file exists and was removed
		"""
		if (sz_file_name == None): return False

		if os.path.exists(sz_file_name) and len(sz_file_name) > 0:
			if os.name == 'nt':
				os.unlink(sz_file_name)
			else:
				cmd = "rm " + sz_file_name
				exist_status = os.system(cmd)
				if (exist_status != 0):
					raise Exception("Fail to remove a file")
			return True
		return False

	def remove_dir(self, path_name):
		if (path_name != None and os.path.isdir(path_name)):
			cmd = "rm -r %s*" % (path_name);
			os.system(cmd)

	def make_path(self, path_name):
		if (not os.path.isdir(path_name) and not os.path.isfile(path_name)):
			cmd = "mkdir -p " + path_name
			os.system(cmd)
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to make a path")

	def is_integer(self, n_value):
		try:
			int(n_value)
			return True
		except ValueError:
			return False

	def is_float(self, n_value):
		try:
			float(n_value)
			return True
		except ValueError:
			return False

	def is_gzip(self, file_name):
		"""
		test if the file name ends in gzip
		"""
		return file_name.endswith(".gz")

	def is_fasta(self, sz_file_name):
		"""
		Test Fata file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: " + sz_file_name))
		b_pass = False
		with open(sz_file_name) as handle:
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp[0] == ">"):
					b_pass = True
					break
				else:
					raise IOError(_("Error: the file is not in FASTA format."))
		if (not b_pass): raise IOError(_("Error: file is not in FASTA format."))

		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError(_("Error: file is not in FASTA format."))

	def read_text_file(self, file_name):
		"""
		read text file and put the result in an vector
		"""
		if (not os.path.exists(file_name)):
			raise IOError(_("Error: file '" + file_name + "' doens't exist."))

		vect_out = []
		with (gzip.open(file_name, 'rt') if file_name.endswith('.gz') else open(file_name)) as handle:
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				vect_out.append(sz_temp)
		return vect_out

	def test_exist_file(self, file_name):
		if (not os.path.exists(file_name)):
			sys.exit("Error: file does not exist - " + file_name)

	def calc_asi(self, rscu_1, rscu_2, vect_codons):
		""" calculation of overall codon usage pattern """
		total_a_b = 0
		square_total_a = 0
		square_total_b = 0
		for codon in vect_codons:
			total_a_b += rscu_1[codon] * rscu_2[codon]
			square_total_a += rscu_1[codon] * rscu_1[codon]
			square_total_b += rscu_2[codon] * rscu_2[codon]
		return ((1.0 - (total_a_b / math.sqrt(square_total_a * square_total_b))) / 2.0)

	def get_aminos(self, b_add_all=False):
		""" get all aminos and codons by amino """
		dt_amino = {}
		for codon in ct.unambiguous_dna_by_id[self.genetic_code].forward_table:
			amino = ct.unambiguous_dna_by_id[self.genetic_code].forward_table[codon]
			if amino in dt_amino:
				dt_amino[amino].append(codon)
			else:
				dt_amino[amino] = [codon]

		if (b_add_all): dt_amino[self.ALL_aminos] = list(ct.unambiguous_dna_by_id[self.genetic_code].forward_table)
		return dt_amino

	def copy_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "cp " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to make a copy a file")

			### set attributes to file 664
			os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

	def move_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "mv " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				raise Exception("Fail to make a copy a file")

	def str2bool(self, v):
		"""
		str to bool
		"""
		return v.lower() in ("yes", "true", "t", "1", "y")

	def get_cai_in_record(self, vect_records, id_, rscu):
		"""  return CAI for an ID in record"""
		for record in vect_records:
			if record.id == id_ and (len(str(record.seq)) % 3) == 0:
				f_cai = CAI(self.clean_fasta_seq(str(record.seq)), RSCUs=rscu)
				if not math.isnan(f_cai): return f_cai
				return None
		return None

	def complement_base(self, base):
		if base == 'A': return 'T'
		if base == 'T': return 'A'
		if base == 'G': return 'C'
		if base == 'C': return 'G'
		return base

	def get_actual_time(self):
		return strftime("%Y-%m-%d %H:%M:%S", localtime())
