from lande.utilities.table import TableFormatter

class PWNFormatter(TableFormatter):
    def pwn(self, pwn):
        pwn = pwn.replace('PSR','')
        if self.table_type == 'latex':
            pwn = pwn.replace('-','$-$')
        elif self.table_type == 'confluence':
            pass
        else:
            raise Exception("...")
        return pwn

