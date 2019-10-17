import random, binascii



class HashOutput(bytearray):
	def hexdigest(self) -> str:
		return binascii.hexlify(self).decode('ascii')


class PearsonHasher:

	def __init__(self, length: int, seed = 'ΑΓΕΩΜΕΤΡΗΤΟΣ ΜΗΔΕΙΣ ΕΙΣΙΤΩ'):
		self.length = length
		generator = random.Random()
		generator.seed(seed)
		self.table = list(range(256))
		generator.shuffle(self.table)


	def hash(self, data: bytes) -> HashOutput:
		result = HashOutput()
		for byte in range(self.length):
			h = self.table[(data[0] + byte) % 256]
			for c in data[1:]:
				h = self.table[h ^ c]
			result.append(h)
		return result



if __name__ == '__main__':
	for L in (2, 4, 10):
		print(L, 'bytes')
		hasher = PearsonHasher(L)
		for s in ('Lorem ipsum dolor sit amet', 'Morem ipsum dolor sit amet', 'Lorem ipsum dolor sit amet!', 'a', 'b'):
			print(s, '\t', hasher.hash(s.encode('utf-8')).hexdigest())
		print()
