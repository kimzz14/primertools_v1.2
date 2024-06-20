class Parameter:
    def __init__(self):
        self.min_primer_GC = 0.4
        self.max_primer_GC = 0.6

        self.min_primer_Tm = 50
        self.max_primer_Tm = 60

        self.min_primer_len = 18
        self.max_primer_len = 25

        self.max_primer_TmDiff = 1.0

        self.min_product_len = 200
        self.max_product_len = 400

        self.is_gc_at_3_prime = True