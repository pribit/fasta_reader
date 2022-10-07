from typing import IO, Optional


def max_gc_by_split():
    max_gc_content = 0
    max_indicator = ''

    with open('max_gc_dataset.txt') as inf:
        fasta = inf.read().split('>')
        del fasta[0]

        for i in fasta:
            fasta_content = i.replace('\n', '')
            indicator = fasta_content[:13]
            dna_string = fasta_content[13:]

            gc_content = (dna_string.count('G') + dna_string.count('C')) / len(dna_string) * 100
            if gc_content > max_gc_content:
                max_gc_content = gc_content
                max_indicator = indicator

        print(max_indicator, max_gc_content)


def max_gc_by_lines():
    max_gc_content = 0
    max_indicator = ''

    with open('max_gc_dataset.txt') as inf:
        for line in inf:
            if line.startswith('>'):
                pass
            else:
                pass


class FastaSequence:
    __slots__ = ('label', 'dna')

    def __init__(self, label: str, dna: str) -> None:
        self.label: str = label
        self.dna: str = dna

    @property
    def gc_content(self) -> float:
        return (self.dna.count('G') + self.dna.count('C')) / len(self.dna) * 100


class FastaReader:
    def __init__(self, fasta_file: IO) -> None:
        self.fasta_file: IO = fasta_file
        self.next_label: Optional[str]
        self.eof: bool

    def __iter__(self) -> 'FastaReader':
        self.fasta_file.seek(0)
        self.next_label = None
        self.eof = False
        return self

    def __next__(self) -> FastaSequence:
        if self.eof:
            raise StopIteration('End of fasta file was reached')
        label: Optional[str] = None or self.next_label
        dna: str = ''

        for line in self.fasta_file:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if not label:
                    label = line.lstrip('>')
                else:
                    self.next_label = line.lstrip('>')
                    return FastaSequence(label=label, dna=dna)
            else:
                dna += line

        self.eof = True
        return FastaSequence(label=label, dna=dna)


def max_gc_by_fasta_reader_iterator():
    max_gc_fasta_sequence: Optional[FastaSequence] = None

    with open('max_gc_dataset.txt') as inf:
        for fasta_sequence in FastaReader(inf):
            if max_gc_fasta_sequence is None or fasta_sequence.gc_content > max_gc_fasta_sequence.gc_content:
                max_gc_fasta_sequence = fasta_sequence

        print(max_gc_fasta_sequence.label, max_gc_fasta_sequence.gc_content)


if __name__ == '__main__':
    max_gc_by_split()
    max_gc_by_fasta_reader_iterator()
