""" Code to clasify regions. """



class PWNClassifier(object):
    """ Classify a PWN """
    def __init__(self, results):
        self.results = results

        at_pulsar_gtlike = results['at_pulsar']['gtlike']
        at_pulsar_pointlike = results['at_pulsar']['pointlike']


        point_gtlike = results['point']['gtlike']
        point_pointlike = results['point']['pointlike']

        extended_gtlike = results['extended']['gtlike']
        extended_pointlike = results['extended']['pointlike']

        self.cutoff=point_gtlike['test_cutoff']

        self.ts_point = max(point_gtlike['TS'],0)
        self.ts_ext = max(extended_gtlike['TS']-point_gtlike['TS'],0)

        try:
            self.ts_cutoff = max(self.cutoff['TS_cutoff'],0)
        except:
            self.ts_cutoff = None

        if self.ts_point > 25:

            if self.ts_cutoff > 16:

                if self.ts_ext > 16:
                    self.type = 'confused'
                    self.spatial = 'at_pulsar'
                else:
                    self.type = 'pulsar'
                    self.spatial = 'point'
            
            else:
                if self.ts_ext > 16:
                    self.type = 'pwn'
                    self.spatial = 'extended'
                else:
                    self.type = 'confused'
                    self.type = 'point'

        else:
            self.type = 'upper_limit'
            self.spatial = 'at_pulsar'
