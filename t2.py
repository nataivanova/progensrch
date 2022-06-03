#!/usr/bin/env python3

import os
import uuid
import logging
import progenlib as p
from telegram.ext import Updater, MessageHandler, Filters

logger = logging.getLogger('progen tool')
logging.basicConfig( level=logging.INFO
                    , format = 'id=' +  req_id + ', t=%(asctime)s, level=%(levelname), message=%(message)s' )

BOT_TOKEN = os.environ.get('TELEGRAM_TOKEN')
db_location = os.environ['DB_LOCATION']

db = p.ProgenitorDatabase(db_location)

def progen_search(req_id):
    query   = p.ProgenitorQuery('/tmp/' + req_id)
    results = p.ProgenitorSearch(query, db)
    s.do_search()
    return s

def downloader(update, context):

    req_id = str(uuid.uuid4())

    with open('/tmp/' + req_id, 'wb') as infile:
        context.bot.get_file(update.message.document).download( out = infile)

    context.bot.send_message(
        chat_id = update.message.chat_id
        , text    = "your search id is " + req_id)
    try:
        s = progen_search(req_id)
        outpath = '/tmp/' + req_id + '.result'
        try:
            outfile = open(outpath, 'wb')
            outfile.write(bytes(str(s),'UTF-8'))
            outfile.close()
            outfile = open(outpath, 'rb')
            context.bot.send_document(
                chat_id  = update.message.chat_id
                , document = outfile
                , disable_content_type_detection = True)
        except Exception as e:
            print (str(e))

    except Exception as e:
        print (e)
        try:
            errpath = '/tmp/' + req_id + '.err'
            errfile = open(errpath, 'wb')
            errfile.write(bytes(str(e),'UTF-8'))
            errfile.close()
            errfile = open(errpath, 'rb')
            context.bot.send_document(
                chat_id = update.message.chat_id
                , document=errfile
                , disable_content_type_detection=True)
        except Exception as e:
            print(str(e))

updater = Updater(BOT_TOKEN, use_context=True)

updater.dispatcher.add_handler(MessageHandler(Filters.document, downloader))

updater.start_polling()
updater.idle()
