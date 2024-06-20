def find_closest_leq(arr, target):
        if not arr:
            return None
        if len(arr) == 1:
            return 0 if arr[0] <= target else None

        low, high = 0, len(arr) - 1

        arrIDX = 0
        while low <= high:
            mid = (low + high) // 2
            if arr[mid] <= target:
                arrIDX = mid
                high = mid - 1
            else:
                low = mid + 1

        return arrIDX

def find_closest_geq(arr, target):
        if not arr:
            return None
        if len(arr) == 1:
            return 0 if arr[0] >= target else None

        low, high = 0, len(arr) - 1

        arrIDX = 0
        while low <= high:
            mid = (low + high) // 2
            if arr[mid] >= target:
                arrIDX = mid
                high = mid - 1
            else:
                low = mid + 1

        return arrIDX


class Primer:
    def __init__(self, seq):
        self.seq = seq.upper()
        self.length = len(self.seq)
        self.gc = 0
        self.tm = 0

        self.hit_DICT = {}

        self.hitSortedFlag_DICT = {}
    
    def add_hit(self, seqName, strand, pos):
        if not seqName in self.hit_DICT: 
            self.hit_DICT[seqName] = {'+':[], '-':[]}
            self.hitSortedFlag_DICT[seqName] = {'+':False, '-':False}

        if strand == '+':
            self.hit_DICT[seqName][strand] += [pos]
        else:
            self.hit_DICT[seqName][strand] += [pos + self.length - 1]
    
    def get_sortedHit_LIST(self, seqName, strand):
        if not seqName in self.hit_DICT:
            return []

        if self.hitSortedFlag_DICT[seqName][strand] == False:
            self.hitSortedFlag_DICT[seqName][strand] = True

            hit_LIST = self.hit_DICT[seqName][strand]
            sortedHit_LIST = sorted(hit_LIST)
            self.hit_DICT[seqName][strand] = sortedHit_LIST

        return self.hit_DICT[seqName][strand]
    
    def get_seqName_SET(self):
        return set(self.hit_DICT.keys())
    
    def text(self):
        return self.seq + ':' + str(self.gc) + ':' + str(self.tm)

class Product:
    def __init__(self, seqName, fTag, rTag, fPrimer, rPrimer, fPos, rPos):
        self.seqName = seqName
        self.fTag = fTag
        self.rTag = rTag
        self.fPrimer = fPrimer
        self.rPrimer = rPrimer
        self.fPos = fPos
        self.rPos = rPos

    def text(self):
        context = []
        context += [self.seqName]
        context += [self.fTag + self.rTag]
        context += [self.fPos]
        context += [self.rPos]
        context += [self.rPos - self.fPos + 1]
        return '\t'.join(map(str, context)) 

class PrimerPair:
    def __init__(self, fPrimer, rPrimer):
        self.fPrimer = fPrimer
        self.rPrimer = rPrimer
        self.product_LIST = []

    def text(self):
        context = []

        if len(self.product_LIST) == 0:
            productN = 'Non-Amplified'
        elif len(self.product_LIST) == 1:
            productN = 'Single-Band'
        else:
            productN = 'Multiple-Bands[' + str(len(self.product_LIST)) + ']'


        common = []
        common += [self.fPrimer.seq]
        common += [self.rPrimer.seq]
        common += [self.fPrimer.length]
        common += [self.rPrimer.length]
        common += ["{:.2f}".format(self.fPrimer.tm)]
        common += ["{:.2f}".format(self.rPrimer.tm)]
        common += ["{:.2f}".format(self.fPrimer.gc)]
        common += ["{:.2f}".format(self.rPrimer.gc)]


        if len(self.product_LIST) == 0:
            return '\t'.join(map(str, common + [productN] + ['NA']*4))

        for product in self.product_LIST:
            context += ['\t'.join(map(str, common + [productN, product.text()]))]

        return '\n'.join(map(str, context))
    
    def find_hit(self, maxSize):
        minSize = self.fPrimer.length + self.rPrimer.length

        seqName_SET = self.fPrimer.get_seqName_SET().intersection(self.rPrimer.get_seqName_SET())
        for seqName in seqName_SET:
            
            for fTag, fPrimer in zip(['F','R'], [self.fPrimer, self.rPrimer]):
                
                fHit_LIST = fPrimer.get_sortedHit_LIST(seqName, '+')
                
                if len(fHit_LIST) == 0 : continue

                for rTag, rPrimer in zip(['F','R'], [self.fPrimer, self.rPrimer]):
                    rHit_LIST = rPrimer.get_sortedHit_LIST(seqName, '-')
                    if len(rHit_LIST) == 0: continue

                    start_fHitIDX = find_closest_leq(fHit_LIST, rHit_LIST[0])
                    if start_fHitIDX == None: continue

                    for fHitIDX in range(start_fHitIDX, len(fHit_LIST)):
                        fPos = fHit_LIST[fHitIDX]

                        start_rHitIDX = find_closest_geq(rHit_LIST, fPos)
                        if start_rHitIDX == None: continue

                        for rHitIDX in range(start_rHitIDX, len(rHit_LIST)):
                            rPos = rHit_LIST[rHitIDX]

                            productSize = rPos - fPos + 1
                            if productSize < minSize: continue
                            if productSize <= maxSize:
                                product = Product(seqName, fTag, rTag, fPrimer, rPrimer, fPos, rPos)
                                self.product_LIST += [product]

                            if maxSize < productSize: break


