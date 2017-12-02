block = """
In our everyday lives we have grown accustomed to computers being a
constant. Computers solve problems that require years for a human in
just seconds. Yet, there is variance in how long it takes a computer
to do the same task multiple times. This variance in performance is a
pervasive problem across computer science. In my talk I will share
three significant ways that performance variance impacts computer
systems in the real world, and introduce you to some modelling methods
that we are using to manage and predict variance in the future.
"""

all_words = {}
to_remove = ',\n'
# Recursive function for removing all characters in 'to_remove'
remove = lambda to_rem, string: (
    remove(to_rem[1:], string.replace(to_rem[0],""))
    if len(to_rem) > 0 else string)

# Correct the entire block with some simple rules
block = block.lower()
block = block.replace("\n"," ")

sentences = []

for sentence in (block.lower()).split("."):
    s_words = {}
    for word in sentence.split(" "):
        word = remove(to_remove, word)
        if len(word) == 0: continue

        all_words[word] = all_words.get(word, 0) + 1
        s_words[word] = s_words.get(word, 0) + 1
    s_words = list(s_words.items())
    s_words.sort(key=lambda k: -k[-1])
    if len(s_words) == 0: continue
    print(s_words)
    print()

all_words = list(all_words.items())
all_words.sort(key=lambda k: -k[-1])
for w in all_words:
    print(w)
