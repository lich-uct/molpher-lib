class FindClosest:

    def __init__(self):
        self.closest = None

    def __call__(self, morph):
        if not self.closest:
            self.closest = morph.copy()
            return
        current_dist = self.closest.getDistToTarget()
        morph_dist = morph.getDistToTarget()
        if morph_dist < current_dist:
            self.closest = morph.copy()