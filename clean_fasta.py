def filter_sequences(input_file, output_file, forbidden_words):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    filtered_sequences = []
    current_sequence = ''
    headers_encountered = set()

    for line in lines:
        if line.startswith('>'):
            header = line.strip()
            keep_sequence = all(word not in header.lower() for word in forbidden_words) and header not in headers_encountered
            if keep_sequence:
                if current_sequence:  # Append the previous sequence if not empty
                    filtered_sequences.append(current_sequence)
                current_sequence = line
                headers_encountered.add(header)
        else:
            if keep_sequence:
                current_sequence += line

    # Append the last sequence after the loop
    if current_sequence and keep_sequence:
        filtered_sequences.append(current_sequence)

    with open(output_file, 'w') as f:
        for sequence in filtered_sequences:
            f.write(sequence)


if __name__ == "__main__":
    input_file = "/path/to/input/fasta"
    output_file = "/path/to/output/fasta"
    forbidden_words = ["like", "putative", "predicted"]
    filter_sequences(input_file, output_file, forbidden_words)
