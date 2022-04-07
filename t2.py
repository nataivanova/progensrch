#!/usr/bin/env python3

import os
import tempfile
import logging

import numpy as np
import os
import pickle

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.DEBUG)

from telegram.ext import Updater, MessageHandler, Filters

BOT_TOKEN = os.environ.get('TELEGRAM_TOKEN')

import uuid
import progenlib as p
import os

def progen_search(req_id):
    x = p.ProgenitorQuery('/tmp/' + req_id)
    db_location = os.environ['DB_LOCATION']
    s = p.ProgenitorSearcher(x.query, db_location)
    s.do_search()
    return s

def downloader(update, context):

    req_id = str(uuid.uuid4())

    with open('/tmp/' + req_id, 'wb') as f:
        context.bot.get_file(update.message.document).download(out=f)
        f.close()
        context.bot.send_message(chat_id=update.message.chat_id, text = "your search id is" + req_id)
        s = progen_search(req_id)
        outpath = '/tmp/' + req_id + '.result'
        try:
            outfile = open(outpath, 'wb')
            outfile.write(bytes(str(s),'UTF-8'))
            outfile.close()
            outfile = open(outpath, 'rb')
            context.bot.send_document(chat_id=update.message.chat_id, document=outfile, disable_content_type_detection=True)
        except Exception as e:
            print (str(e))


    # tempfile.NamedTemporaryFile(
          # mode   = 'r+t'
        # , suffix = '.txt'
        # , prefix = 'progsrch.up.'
        # , delete = False
        # , errors = None)  as f:
        # context.bot.get_file(update.message.document).download(out=f)
        # print("downloaded tp file" + f + "\n")

updater = Updater(BOT_TOKEN, use_context=True)

updater.dispatcher.add_handler(MessageHandler(Filters.document, downloader))

updater.start_polling()
updater.idle()
