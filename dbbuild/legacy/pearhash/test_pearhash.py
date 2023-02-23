import unittest
from pearhash import PearsonHasher



class TestPearsonHasher(unittest.TestCase):

	def test_table_is_a_permutation_of_range_256(self):
		hasher = PearsonHasher(2)
		self.assertEqual(set(hasher.table), set(range(256)))


	def test_two_bytes(self):
		hasher = PearsonHasher(2)
		self.assertEqual(hasher.hash(b'ni hao').hexdigest(), '1297')


	def test_two_bytes_custom_seed(self):
		hasher = PearsonHasher(2, seed = 'whatevs')
		self.assertEqual(hasher.hash(b'ni hao').hexdigest(), 'd710')


	def test_four_bytes(self):
		hasher = PearsonHasher(4)
		self.assertEqual(hasher.hash(b'ni hao').hexdigest(), '1297b8d9')
