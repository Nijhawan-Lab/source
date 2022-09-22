from Bio import Seq
import jpbio.util
import matplotlib.pyplot as plt

class QualityStatistics:

    def __init__(self):
        self.read_counts = []
        self.read_lengths = []
        self.quality_scores = []
        self.N_counts = []
        self.N_sequences = 0

    def addSequence(self, read):
        quality = read.letter_annotations["phred_quality"]
        seq_length = len(read)
        self.read_counts = jpbio.util.pad_list(self.read_counts, seq_length, 0)
        self.N_counts = jpbio.util.pad_list(self.N_counts, seq_length, 0)
        self.quality_scores = jpbio.pad_dict_list(self.quality_scores, seq_length)
        self.read_lengths.append(seq_length)

        for i in range(seq_length):
            self.read_counts[i] += 1
            if read.seq[i] == 'N':
                self.N_counts[i] += 1
            q = read.letter_annotations["phred_quality"][i]
            if q in self.qualtiy_scores[i]:
                self.qualtiy_scores[i][ q ] += 1
            else:
                self.qualtiy_scores[i][ q ] = 1
        
        self.N_sequences += 1
        return self.N_sequences

    def DiagnosticPlot(self, maintitle="Read Quality"):

        def mean_from_dictionary_of_counts(d):
            M = 0
            T = 0
            for k in d:
                M += k*d[k]
                T += d[k]
            return M/T
        
        max_length = len(self.read_counts)
        positions = range(1, max_length + 1)
        max_qual_by_position = [max( [k for k in d] ) for d in self.quality_scores]
        min_qual_by_position = [min( [k for k in d] ) for d in self.quality_scores]
        mean_qual_by_position = [mean_from_dictionary_of_counts(d) for d in self.quality_scores]
        max_N_count = max([10, max( self.N_counts )])
        
        fig, (axis_quality, axis_N, axis_read_count) = plt.subplots(
            3,1,
            sharex='all',
            figsize=(15,7)
        )
        fig.suptitle(maintitle, fontsize=16)
        
        min_max_handle = axis_quality.fill_between(
            x=positions, 
            y1=min_qual_by_position, 
            y2=max_qual_by_position, 
            color="lightgrey", 
            label='range')
        mean_handle, = axis_quality.step(
            x=positions,
            y=mean_qual_by_position, 
            label='mean')
        axis_quality.legend(loc='lower center')
        axis_quality.set_xlabel('position')
        axis_quality.set_ylabel('quality')
        axis_quality.yaxis.grid(True)
        axis_quality.spines["bottom"].set_visible(True)
        axis_quality.spines["top"].set_visible(False)
        axis_quality.spines["left"].set_visible(True)
        axis_quality.spines["right"].set_visible(False)

        axis_read_count.set_yscale('log')
        axis_read_count.plot(self.read_counts, color='grey')
        axis_read_count.hist(self.read_counts, range=(0,max_length), bins=25, color='lightblue')
        axis_read_count.set_ylabel('read count')
        axis_read_count.spines["bottom"].set_visible(True)
        axis_read_count.spines["top"].set_visible(False)
        axis_read_count.spines["left"].set_visible(True)
        axis_read_count.spines["right"].set_visible(False)
        
        axis_N.spines["bottom"].set_visible(True)
        axis_N.spines["top"].set_visible(False)
        axis_N.spines["left"].set_visible(True)
        axis_N.spines["right"].set_visible(False)
        axis_N.set_ylim([0, max_N_count])
        axis_N.step(x=positions, y=self.N_counts)
        axis_N.set_ylabel('N content')

        return fig
