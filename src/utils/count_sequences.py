"""Count of excluded and included sequences"""


class CountSequences:

    def __init__(self):
        self.count_divisible_3 = 0  # sequences that are not divisible per 3
        self.count_pass = 0         # sequences that pass to the teste

    def add_divisible_3(self):
        self.count_divisible_3 += 1

    def add_pass(self):
        self.count_pass += 1

    def __str__(self):
        return f"Fail Divisible 3: {self.count_divisible_3}\nPassed genes: {self.count_pass}"
