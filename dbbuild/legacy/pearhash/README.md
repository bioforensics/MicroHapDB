pearhash
========

Simple [Pearson hashing](http://en.wikipedia.org/wiki/Pearson_hashing) algorithm. Python 3 only.


Instalation
-----------

```
pip install pearhash
```


Usage
-----

```
>>> from pearhash import PearsonHasher
>>> hasher = PearsonHasher(2) # Set desired hash length in bytes.
>>> hasher.hash(b'ni hao')
bytearray(b'\x12\x97')
>>> hasher.hash(b'ni hao').hexdigest()
'1297'
```

Hash length can be changed easily:

```
>>> hasher = PearsonHasher(4)
>>> hasher.hash(b'ni hao').hexdigest()
'1297b8d9'
```

And you can even change the seed used to generate the pseudorandom table:

```
>>> hasher = PearsonHasher(2, seed = 'whatevs')
>>> hasher.hash(b'ni hao').hexdigest()
'd710'
```
