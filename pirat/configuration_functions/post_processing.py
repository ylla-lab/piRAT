import random
import os

ALL_QUOTES = [
    # Bioinformatics & Science Quotes
    "Bioinformatics is like assembling IKEA furniture without instructions or all the screws.\n\t- Dr. A. Witty, Genomics Conference 2019",
    "Sequencing teaches us that what you get is not always what you expected.\n\t- Adapted from J. Craig Venter",
    "Every dataset is a new puzzle, and sometimes you don’t even know what the picture is.\n\t- Anonymous PI, Lab Humor Slack, 2021",
    "The best part of bioinformatics is making sense out of biological chaos.\n\t- Dr. Marie Data, 'Reflections on Science', 2020",
    "Good code is its own documentation, but bad code is a lifelong mystery.\n\t- Based on Phil Karlton",
    "The genome is the world’s longest palindrome, if you squint hard enough.\n\t- Internet meme, 2022",
    "The only thing that moves faster than technological advances is the rate of new sequencing errors.\n\t- Dr. Sequence R. Fast, ISMB 2022",
    "If at first you don’t succeed, blame batch effects.\n\t- Lab Wisdom, University of Nowhere",
    "NGS: Because we wanted bigger problems to solve.\n\t- Dr. Ada Byte, 'Computing Tomorrow', 2021",
    "In science, you never finish—your funding just runs out.\n\t- Common grad student saying",

    # Programming & Tech Quotes
    "To code is human, to debug is divine.\n\t- Frank Savage, paraphrased",
    "The most effective debugging tool is still careful thought, coupled with judiciously placed print statements.\n\t- Brian Kernighan, 'The Elements of Programming Style', 1978",
    "There are two ways to write error-free programs; only the third one works.\n\t- Alan J. Perlis",
    "If at first you don’t succeed, call it version 1.0.\n\t- Programmer’s proverb",
    "Programming is thinking, not typing.\n\t- Casey Patton",
    "Any fool can write code that a computer can understand. Good programmers write code that humans can understand.\n\t- Martin Fowler, 'Refactoring', 1999",
    "Weeks of coding can save you hours of planning.\n\t- Anonymous",
    "Software and cathedrals are much the same: first we build them, then we pray.\n\t- Sam Redwine, 1995",
    "There is no place like 127.0.0.1.\n\t- T-shirt wisdom",
    "If debugging is the process of removing software bugs, then programming must be the process of putting them in.\n\t- Edsger W. Dijkstra",

    # Science & Life Quotes
    "Research is what I’m doing when I don’t know what I’m doing.\n\t- Wernher von Braun",
    "Science is the belief in the ignorance of experts.\n\t- Richard Feynman, 1966",
    "If we knew what we were doing, it wouldn’t be called research.\n\t- Albert Einstein",
    "The good thing about science is that it’s true whether or not you believe in it.\n\t- Neil deGrasse Tyson",
    "It’s not that I’m so smart, it’s just that I stay with problems longer.\n\t- Albert Einstein",
    "In theory, theory and practice are the same. In practice, they are not.\n\t- Yogi Berra",
    "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.\n\t- Stephen Hawking",
    "Life is like riding a bicycle. To keep your balance, you must keep moving.\n\t- Albert Einstein",

    # Pop Culture & Movies
    "I have no special talents. I am only passionately curious.\n\t- Albert Einstein",
    "Do, or do not. There is no try.\n\t- Yoda, Star Wars: The Empire Strikes Back (1980)",
    "The truth is out there.\n\t- The X-Files (TV Series, 1993)",
    "It’s a trap!\n\t- Admiral Ackbar, Star Wars: Return of the Jedi (1983)",
    "Just keep swimming.\n\t- Dory, Finding Nemo (2003)",
    "Houston, we have a problem.\n\t- Jim Lovell, Apollo 13 (1995)",
    "I’m not superstitious, but I am a little stitious.\n\t- Michael Scott, The Office (TV Series, 2007)",
    "Elementary, my dear Watson.\n\t- Sherlock Holmes (paraphrase, not in Conan Doyle books)",
    "The cake is a lie.\n\t- Portal (Video Game, 2007)",
    "Winter is coming.\n\t- Game of Thrones (TV Series, 2011)"
]



def display_random_quote():
    """Displays a random quote from the combined list"""
    quote = random.choice(ALL_QUOTES)
    print(f"\n\"{quote}\"")

