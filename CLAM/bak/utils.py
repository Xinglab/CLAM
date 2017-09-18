
class CLAM_mread(object):
	"""
	object to store a read alignment
	"""
	def __init__(self, alignment, read_tagger_func, tag_of_interest = ['NH']):
		self.reference_id = alignment.reference_id
		self.is_reverse = alignment.is_reverse
		self.cigarstring = alignment.cigarstring
		self.pos = alignment.pos
		self.qname = alignment.qname
		self.tag = read_tagger_func(alignment)
		self.flag = alignment.flag
		self.alignment_tags = [x for x in alignment.tags if x[0] in tag_of_interest]
	
	def __eq__(self, other):
		this = (self.reference_id, self.is_reverse, self.pos, self.qname)
		that = (other.reference_id, other.is_reverse, other.pos, other.qname)
		return this == that
	
	def __hash__(self):
		return hash((self.reference_id, self.is_reverse, self.pos, self.qname))
	
	def __str__(self):
		s = "\t".join([
			self.qname, 
			str(self.flag),
			str(self.reference_id), 
			str(self.pos),
			self.cigarstring])
		return s