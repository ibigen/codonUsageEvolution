
import os
import socket

class Expression(object):
    def __init__(self):
        file_information = "E.coli_information.txt"
        file_expression = "E.coli_expression.txt"
        if (socket.gethostname() == "cs-nb0008"):  # test computer name
            self.base_path = "/home/projects/ua/master/codon_usage"
        else:
            self.base_path = r"C:\Users\Francisca\Desktop\TeseDeMestrado"

        self.file_information = os.path.join(self.base_path, file_information)
        self.file_expression = os.path.join(self.base_path, file_expression)


    def samples_information(self):
        information = []
        with open(self.file_information, 'r') as information_file:
            for line in information_file:
                line = line.strip()
                line = line.split('\t')
                information.append(line)


    def expression_values(self):
        expression = []
        with open(self.file_expression, 'r') as expression_file:
            for line in expression_file:
                line = line.strip()
                line = line.split('\t')
                expression.append(line)



